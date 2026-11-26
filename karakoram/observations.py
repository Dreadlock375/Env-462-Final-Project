"""Prepare observed daily discharge independently of model results."""
from __future__ import annotations
import collections, csv, datetime as dt, math, re
from pathlib import Path
from .common import ROOT, OBS, CONFIG, require, sha, read_json, atomic_json, write_csv
STATION_ID = 'Indus_at_Besham_Qila'
START = dt.date(2005, 1, 1)
END = dt.date(2024, 12, 31)
CLEANED_LOCATIONS = (
    Path('data/observations/wapda/original/WAPDA/Cleaned_Data'),
    Path('data/observations/wapda/cleaned'),
)
OUTPUT_RELATIVE = Path('data/observations/wapda/processed')
DEFAULT_PERIODS = {'calibration': (2005, 2012),
                   'validation': (2013, 2019),
                   'holdout': (2020, 2024)}
DATE_HEADERS = ('date', 'day', 'datetime', 'timestamp', 'observation_date')
VALUE_HEADERS = ('discharge_m3s', 'discharge_m3_s', 'discharge_m3_per_s',
                 'flow_m3s', 'flow_m3_s', 'flow_m3_per_s', 'q_m3s',
                 'observed_discharge_m3s', 'discharge_cumecs', 'flow_cumecs',
                 'discharge', 'flow', 'value_m3s', 'value_m3_s', 'obs_m3s',
                 'discharge_m3sec', 'value', 'm3s')
FLAG_HEADERS = ('flag', 'qc_flag', 'quality_flag', 'source_flag', 'measurement_flag')
MISSING_WORDS = {'', '-', '--', 'na', 'n/a', 'nan', 'null', 'none', 'missing',
                 'data missing', '#ref!', '-1', '-1.000', '-9999'}
EXPECTED_STATION_CODE = '34729801'

def normalized_header(header: str) -> str:
    """Normalize minor separators and unit notation without guessing semantics."""
    s = header.strip().lower().replace('³', '3').replace(' ', '_')
    s = s.replace('(', '_').replace(')', '').replace('/', '_')
    return re.sub(r'_+', '_', re.sub(r'[^a-z0-9_]+', '_', s)).strip('_')

def resolve_column(headers: list[str], aliases: tuple[str, ...], label: str) -> str | None:
    """Resolve one column; fail rather than silently choose ambiguous columns."""
    available = collections.defaultdict(list)
    for h in headers:
        available[normalized_header(h)].append(h)
    for alias in aliases:
        found = available.get(normalized_header(alias), [])
        if len(found) > 1:
            raise ValueError(f'Ambiguous {label} headers: {found}')
        if found:
            return found[0]
    return None

def read_date(value: str) -> dt.date:
    """Parse unambiguous cleaned timestamps, rejecting date/time shifts."""
    s = value.strip()
    for fmt in ('%Y-%m-%d', '%Y/%m/%d', '%Y-%m-%d %H:%M:%S',
                '%Y-%m-%dT%H:%M:%S', '%d-%b-%Y', '%d/%m/%Y'):
        try:
            timestamp = dt.datetime.strptime(s, fmt)
            if timestamp.time() != dt.time():
                raise ValueError(f'Non-midnight timestamp not accepted as a daily observation: {s}')
            return timestamp.date()
        except ValueError as exc:
            if 'Non-midnight' in str(exc):
                raise
    raise ValueError(f'Unsupported or ambiguous date {s!r}; cleaned ISO dates are preferred')

def read_discharge(raw: str) -> float | None:
    """Interpret a cleaned m³/s value; never interpolate or coerce bad measurements."""
    s = raw.strip()
    if s.lower() in MISSING_WORDS:
        return None
    if ',' in s:
        if not re.fullmatch(r'[+-]?\d{1,3}(,\d{3})+(\.\d+)?', s):
            raise ValueError(f'Ambiguous comma-formatted discharge {s!r}')
        s = s.replace(',', '')
    try:
        value = float(s)
    except ValueError as exc:
        raise ValueError(f'Unrecognized nonnumeric discharge {s!r}') from exc
    if not math.isfinite(value) or value < 0:
        raise ValueError(f'Nonfinite or negative discharge {s!r}')
    return value

def read_daily_csv(path: Path, column_overrides: dict[str, str] | None = None) -> tuple[dict[dt.date, dict], dict]:
    """Parse source; return records plus source-header and QA metadata."""
    with path.open('r', newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames
        if not headers or len(headers) != len(set(headers)):
            raise ValueError('Missing or duplicated CSV headers')
        column_overrides = column_overrides or {}
        for kind, chosen in column_overrides.items():
            if chosen not in headers:
                raise ValueError(f'Explicit {kind} column {chosen!r} is not in headers: {headers}')
        date_col = column_overrides.get('date') or resolve_column(headers, DATE_HEADERS, 'date')
        value_col = column_overrides.get('value') or resolve_column(headers, VALUE_HEADERS, 'discharge')
        flag_col = column_overrides.get('flag') or resolve_column(headers, FLAG_HEADERS, 'flag')
        if not date_col or not value_col:
            raise ValueError('Could not identify date/discharge columns. Available headers: ' + repr(headers))
        if not flag_col:
            raise ValueError('Per-station SWHP CSV has no flag column; use the flagged per-station source')
        unit_cols = [h for h in headers if normalized_header(h) in ('unit', 'units')]
        records: dict[dt.date, dict] = {}
        ignored_outside = 0
        first, last = None, None
        missing_in_source = 0
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f'Malformed CSV: extra columns on line {line_number}')
            if all(not (v or '').strip() for v in row.values()):
                continue
            try:
                date = read_date(row[date_col] or '')
            except ValueError as exc:
                raise ValueError(f'Line {line_number}: {exc}') from exc
            # Count excluded dates, but never silently import observations from other years.
            if not START <= date <= END:
                ignored_outside += 1
                continue
            if date in records:
                raise ValueError(f'Duplicate source observation date {date} at line {line_number}')
            if unit_cols:
                units = (row[unit_cols[0]] or '').strip().lower().replace('³', '3').replace(' ', '')
                if units and units not in ('m3/s', 'm3s-1', 'm^3/s', 'cumec', 'cumecs'):
                    raise ValueError(f'Line {line_number}: unexpected discharge unit {units!r}')
            try:
                value = read_discharge(row[value_col] or '')
            except ValueError as exc:
                raise ValueError(f'Line {line_number}: {exc}') from exc
            flag = (row[flag_col] or '').strip()
            if value is None:
                missing_in_source += 1
            elif flag.lower() in {'invalid', 'excluded', 'zero-as-missing', 'missing-code', 'implausible', 'blank', 'month-missing'}:
                raise ValueError(f'Line {line_number}: numerical discharge with invalid/missing flag {flag!r}')
            records[date] = {'date': date.isoformat(), 'discharge_m3s': value,
                             'source_flag': flag, 'measurement_flag': flag.upper() == 'M',
                             'source_line': line_number,
                             'source_fields': {h: row[h] or '' for h in headers}}
            first = date if first is None else min(first, date)
            last = date if last is None else max(last, date)
    if not records:
        raise ValueError('No valid dated records from 2005–2024 found in source')
    meta = {'source_columns': headers, 'date_column': date_col,
            'discharge_column': value_col, 'flag_column': flag_col,
            'n_source_records_in_period': len(records),
            'n_source_missing_values': missing_in_source,
            'n_source_rows_outside_period': ignored_outside,
            'source_first_date': first.isoformat(), 'source_last_date': last.isoformat(),
            'unit_assumption': 'SWHP daily discharge in m3/s per WAPDA Cleaned_Data/DATA_SUMMARY.md'}
    return records, meta


SOURCE = ROOT / 'data/observations/wapda/original/WAPDA/Cleaned_Data/SWHP_discharge/Indus_at_Besham_Qila_daily_discharge.csv'

def prepare():
    records, meta=read_daily_csv(SOURCE)
    expected=[START+dt.timedelta(days=i) for i in range((END-START).days+1)]
    require(sorted(records)==expected and all(records[d]['discharge_m3s'] is not None for d in expected),'Observation coverage incomplete')
    rows=[{'date':d.isoformat(),'discharge_m3s':records[d]['discharge_m3s'],'source_flag':records[d]['source_flag']} for d in expected]
    require(sha(SOURCE)==CONFIG['observation_source_sha256'],'WAPDA source differs from frozen input')
    OBS.mkdir(parents=True,exist_ok=True)
    write_csv(OBS/'daily.csv',rows)
    # Station identity and units recorded in the original preparation report;
    # no inference from coordinates or changing source flags.
    atomic_json(OBS/'report.json',dict(source_sha256=sha(SOURCE),daily_sha256=sha(OBS/'daily.csv'),station=CONFIG['observation_station'],metadata=meta,days=len(rows),limitations=['Gauge equivalence remains provisional.']))
    print('Prepared observations:',OBS)

if __name__=='__main__':prepare()
