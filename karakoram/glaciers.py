"""Acquire calibrated ERA5 glacier directories and produce resumable daily water series."""
from __future__ import annotations
import argparse, calendar, csv, fcntl, json, math, os, re, shutil, tempfile, traceback
from pathlib import Path
from datetime import date, timedelta, datetime, timezone
import numpy as np
from .common import ROOT, CONFIG, BASIN, GLACIERS, require, sha, atomic_json, START, END, DAYS
from . import glacier_physics as physics
from .basin import read_jobs
CLIMATE=Path('data/climate/era5/processed/daily_2005_2019')
NATIVE=Path('native_025deg/era5_daily_2005_2019.nc')
WORK=Path('working/glaciers')
OUT=GLACIERS.relative_to(ROOT)
URL=CONFIG['era5_archive']
INPUTS=('mb_calib.json','climate_historical.nc','model_flowlines.pkl')
METHOD='ERA5_MONTHLY_CONSTRAINED_FIXED_GEOMETRY_V1'
EXPECTED=CONFIG['expected_glaciers']
MIN_FREE=8*1024**3
STATUS='ERA5_GLACIER_DAILY_PREPARED_NOT_INDEPENDENTLY_VALIDATED'
DATES=[d.isoformat() for d in DAYS]
MONTHS=[f'{y}-{m:02}' for y in range(2005,2020) for m in range(1,13)]

def exact_index(coords,value):
    matches=np.flatnonzero(np.isclose(coords,value,rtol=0,atol=1e-5))
    require(len(matches)==1,f'OGGM reference coordinate {value} is not uniquely available on native ERA5 grid')
    return int(matches[0])

def finite(a):
    require(not np.ma.getmaskarray(a).any(), 'Masked climate values')
    a = np.asarray(a, dtype=float)
    require(np.isfinite(a).all(), 'Nonfinite climate values')
    return a

def dates(t):
    from netCDF4 import num2date
    return [date(x.year,x.month,x.day) for x in num2date(t[:],t.units,calendar=getattr(t,'calendar','standard'))]

def load_climate():
    from netCDF4 import Dataset
    base = ROOT / CLIMATE
    p = base / NATIVE
    r = json.loads((base / 'conversion_report.json').read_text())
    digest = sha(p)
    require(r['days'] == len(DAYS) and r['output_sha256'][str(NATIVE)] == digest, 'ERA5 native file hash/coverage mismatch')
    with Dataset(p) as ds:
        require(dates(ds['time']) == DAYS, 'ERA5 dates differ from expected 2005–2019')
        lat = finite(ds['lat'][:])
        lon = finite(ds['lon'][:])
        require(np.array_equal(lat, np.arange(38, 29.99, -0.25)) and np.array_equal(lon, np.arange(71.5, 83.01, 0.25)), 'Unexpected native ERA5 grid')
        fields = {}
        for v, u in [('tas', 'K'), ('pr', 'kg m-2 s-1')]:
            require(ds[v].dimensions == ('time', 'lat', 'lon') and ds[v].units == u, f'Invalid {v} schema')
            fields[v] = finite(ds[v][:])
        require((fields['pr'] >= 0).all(), 'Negative ERA5 precipitation')
    return (lat, lon, fields, digest)

def reference(path, climate):
    from netCDF4 import Dataset
    lat, lon, grid, _ = climate
    with Dataset(path) as ds:
        source = str(ds.getncattr('climate_source'))
        require(source.strip().upper() == 'ERA5', f'Unexpected monthly climate source: {source}')
        y = float(ds.getncattr('ref_pix_lat'))
        x = float(ds.getncattr('ref_pix_lon')) % 360
        i = exact_index(lat, y)
        j = exact_index(lon, x)
        monthly_dates = dates(ds['time'])
        t = finite(ds['temp'][:])
        p = finite(ds['prcp'][:])
        h = float(ds.getncattr('ref_hgt'))
    require(np.isfinite(h), 'Invalid reference elevation')
    raw_t = grid['tas'][:, i, j] - 273.15
    raw_p = grid['pr'][:, i, j] * 86400
    by_month = {}
    monthly = {}
    for d, tt, pp in zip(DAYS, raw_t, raw_p):
        by_month.setdefault((d.year, d.month), []).append((d, float(tt), float(pp)))
    for d, tt, pp in zip(monthly_dates, t, p):
        k = (d.year, d.month)
        if k in by_month:
            require(k not in monthly, 'Duplicated OGGM monthly climate record')
            monthly[k] = (float(tt), float(pp))
    require(set(monthly) == set(by_month), 'Missing OGGM monthly climate record')
    comparison = []
    for k, rows in by_month.items():
        rt = float(np.mean([r[1] for r in rows]))
        rp = float(np.sum([r[2] for r in rows]))
        mt, mp = monthly[k]
        comparison.append({'month': f'{k[0]}-{k[1]:02d}', 'raw_daily_mean_C': rt, 'oggm_uncalibrated_monthly_C': mt, 'temperature_difference_C': rt - mt, 'raw_daily_sum_mm': rp, 'oggm_uncalibrated_monthly_mm': mp, 'precipitation_difference_mm': rp - mp})
    return (by_month, {'source': source, 'latitude': y, 'longitude': x, 'reference_elevation_m': h, 'monthly_raw_comparison': comparison, 'max_abs_raw_temperature_difference_C': max((abs(x['temperature_difference_C']) for x in comparison)), 'max_abs_raw_precipitation_difference_mm': max((abs(x['precipitation_difference_mm']) for x in comparison))})

def utc(): return datetime.now(timezone.utc).isoformat()

def climate_gate(ref):
    # Screening tolerances, not a calibration or statistical validation.
    require(ref['source'].strip().upper()=='ERA5','Non-ERA5 reference')
    rows=ref['monthly_raw_comparison']
    require([r['month'] for r in rows]==MONTHS,'Missing/duplicate monthly reference')
    for r in rows:
        t=r['temperature_difference_C']; d=r['precipitation_difference_mm']; q=r['oggm_uncalibrated_monthly_mm']
        require(all(math.isfinite(x) for x in (t,d,q)) and q>=0,'Invalid monthly reference')
        require(abs(t)<=0.05 and abs(d)<=max(0.5,0.005*q),f'Raw ERA5 monthly mismatch {r["month"]}: {t} C, {d} mm')

def verify_report(r, job, provenance):
    require(r['status'] == STATUS and r['rgi_id'] == job['rgi_id'] and (r['provenance'] == provenance), 'Output provenance mismatch')
    require(r['worklist_entry'] == job and r['daily_days'] == len(DATES), 'Worklist or date mismatch')
    require(r['calibration']['baseline_climate_source'] == 'ERA5', 'Non-ERA5 calibration')
    climate_gate(r['reference_climate'])
    s = r['diagnostic_summary']
    require(s['daily_days'] == len(DATES), 'Invalid diagnostic day count')
    require(s['monthly_absolute_tolerance_m3'] == 1e-06, 'Unexpected closure tolerance')
    for k in ('total_melt_m3', 'rain_m3', 'snowfall_m3'):
        require(0 <= s['monthly_max_abs_residual_m3'][k] <= 1e-06, f'Closure failed: {k}')
    require(math.isfinite(r['glacier_area_m2']) and r['glacier_area_m2'] > 0, 'Invalid modeled area')
    require([x['month'] for x in r['monthly_diagnostics']] == MONTHS, 'Incomplete monthly diagnostics')
    for k in ('total_reference_melt_m3', 'total_reference_rain_m3', 'total_disaggregated_melt_m3', 'total_disaggregated_rain_m3'):
        require(math.isfinite(s[k]) and s[k] >= 0, f'Invalid total: {k}')
    return s

def verify_existing(out, job, prov):
    r = json.loads((out / 'report.json').read_text())
    verify_report(r, job, prov)
    require(sha(out / 'daily.csv') == r['daily_csv_sha256'], 'Daily CSV checksum mismatch')
    return r

def process(job, climate, prov, workflow, utils, model_class):
    gid = job['rgi_id']
    w = ROOT / WORK
    out = ROOT / OUT / gid
    source = w / 'per_glacier' / gid[:8] / gid[:11] / gid
    if out.exists():
        return ('REUSED', verify_existing(out, job, prov), source)
    if not all(((source / f).is_file() for f in INPUTS)):
        obtained = workflow.init_glacier_directories([gid], from_prepro_level=3, prepro_base_url=URL, prepro_border=80, prepro_rgi_version='62')
        require(len(obtained) == 1 and Path(obtained[0].dir).resolve() == source.resolve(), 'Unexpected ERA5 directory location')
    require(all(((source / f).is_file() for f in INPUTS)), 'Incomplete ERA5 glacier directory')
    hashes = {f: sha(source / f) for f in INPUTS}
    cal = json.loads((source / 'mb_calib.json').read_text())
    require(cal['baseline_climate_source'] == 'ERA5' and cal['rgi_id'] == gid, 'Wrong calibration source or glacier ID')
    by_month, ref = reference(source / 'climate_historical.nc', climate)
    climate_gate(ref)
    gd = utils.GlacierDirectory(gid, base_dir=str(w / 'per_glacier'))
    require(Path(gd.dir).resolve() == source.resolve(), 'GlacierDirectory path mismatch')
    heights, areas = physics.glacier_bin_geometry(gd.read_pickle('model_flowlines'))
    model = model_class(gd)
    fields = []
    for y in range(2005, 2020):
        require(model.is_year_valid(y), f'No valid monthly OGGM climate for {y}')
        for m in range(1, 13):
            fields.append(physics.climate_month(model, y, m, by_month[y, m], heights, utils.date_to_floatyear))
    daily, monthly, stats = physics.glacier_daily(fields, areas)
    require([x['date'] for x in daily] == DATES, 'Output date mismatch')
    require(hashes == {f: sha(source / f) for f in INPUTS}, 'OGGM source changed during computation')
    r = {'status': STATUS, 'rgi_id': gid, 'provenance': prov, 'production_runner_sha256': sha(Path(__file__)), 'model_source_sha256': hashes, 'reference_climate': ref, 'calibration': cal, 'glacier_area_m2': float(areas.sum()), 'elevation_bins': len(areas), 'diagnostic_summary': stats, 'monthly_diagnostics': monthly, 'worklist_entry': job, 'daily_days': len(daily), 'warnings': ['Fixed geometry, monthly melt and precipitation constraints imposed, not independent validation.', 'Zero initial snow/firn in 2005; no refreezing, retention, glacier evolution or routing.', 'Daily peaks and basin membership remain unverified.']}
    verify_report(r, job, prov)
    stage = Path(tempfile.mkdtemp(prefix='.' + gid + '-', dir=out.parent))
    try:
        with (stage / 'daily.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=physics.DAILY_FIELDS)
            writer.writeheader()
            writer.writerows(daily)
        r['daily_csv_sha256'] = sha(stage / 'daily.csv')
        (stage / 'report.json').write_text(json.dumps(r, indent=2, allow_nan=False, default=lambda o: o.item() if isinstance(o, np.generic) else str(o)) + '\n')
        require(not out.exists(), 'Output appeared during run; refusing overwrite')
        os.rename(stage, out)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return ('COMPUTED', r, source)

def prune(source, job):
    gid = job['rgi_id']
    base = ROOT / WORK / 'per_glacier'
    require(source == base / gid[:8] / gid[:11] / gid, 'Unexpected prune target')
    for node in [source, *source.parents]:
        require(not node.is_symlink(), f'Refusing symlinked prune target: {node}')
        if node == ROOT:
            break
    if source.exists():
        shutil.rmtree(source)

def summary(entries, total, prov, finished):
    ok = [x for x in entries if x['status'] in ('COMPUTED', 'REUSED')]
    bad = [x for x in entries if x['status'] == 'BLOCKED']
    return {'status': 'ERA5_BASIN_DAILY_PRODUCTION_COMPLETE_PROVISIONAL' if finished and len(ok) == EXPECTED else 'ERA5_BASIN_DAILY_PRODUCTION_INCOMPLETE', 'updated_utc': utc(), 'provenance': prov, 'production_runner_sha256': sha(Path(__file__)), 'expected_glaciers': EXPECTED, 'selected_glaciers': total, 'completed': len(ok), 'computed': sum((x['status'] == 'COMPUTED' for x in ok)), 'reused': sum((x['status'] == 'REUSED' for x in ok)), 'blocked': len(bad), 'not_processed_in_selection': total - len(entries), 'rgi_area_km2_completed': math.fsum((x['rgi_area_km2'] for x in ok)), 'modeled_area_km2_completed': math.fsum((x['modeled_area_km2'] for x in ok)), 'total_melt_m3': math.fsum((x['melt_m3'] for x in ok)), 'total_rain_m3': math.fsum((x['rain_m3'] for x in ok)), 'blocked_entries': bad, 'entries': entries, 'warnings': ['Not independently validated. All release locations and basin membership remain provisional.', 'Totals refer only to successful glaciers in this run; no missing-glacier infilling.']}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fetch-prepro', action='store_true', help='Permit downloading official ERA5 L3 directories')
    ap.add_argument('--prune-gdirs', action='store_true', help='Remove extracted ERA5 directories only after output verification; never remove archives')
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--reproduce-omissions', action='store_true', help='Explicitly omit the seven glaciers excluded from the original validated run')
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    require(a.fetch_prepro, 'Explicit --fetch-prepro required for basin production')
    base = ROOT / OUT
    base.mkdir(parents=True, exist_ok=True)
    with (base / '.production.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        run(a, base)

def run(a,base):
    jobs=read_jobs(ROOT)
    require(len(jobs)==EXPECTED and all(re.fullmatch(r'RGI60-(13|14)\.\d{5}',j['rgi_id']) for j in jobs),'Unexpected worklist')
    require(all(math.isfinite(float(j['rgi_area_km2'])) and float(j['rgi_area_km2'])>0 for j in jobs),'Invalid inventory area')
    # Group by the OGGM archive prefix to maximize cache reuse.
    jobs=sorted(jobs,key=lambda j:j['rgi_id'])
    climate=load_climate(); prov=provenance(climate[3],a.reproduce_omissions)
    import oggm
    from oggm import cfg,workflow,utils
    from oggm.core.massbalance import MonthlyTIModel
    require(oggm.__version__=='1.6.3',f'Unexpected OGGM {oggm.__version__}')
    cfg.initialize_minimal(logging_level='WARNING'); cfg.PARAMS['border']=80;cfg.PARAMS['use_multiprocessing']=False;cfg.PARAMS['baseline_climate']='ERA5';cfg.PARAMS['rgi_version']='62'
    w=ROOT/WORK;w.mkdir(parents=True,exist_ok=True);cfg.PATHS['working_dir']=str(w)
    entries=[]; consecutive=0; finished=False; report=base/'production_report.json';errors=base/'production_errors.log'
    print(f'ERA5 production: {len(jobs)} glaciers; climate loaded once',flush=True)
    try:
        for i,j in enumerate(jobs,1):
            gid=j['rgi_id']; out=base/gid
            if a.reproduce_omissions and gid in CONFIG['reference_omissions']:
                entries.append({'rgi_id':gid,'status':'BLOCKED','error_class':'ReferenceOmission','error':'Explicit historical reproduction omission'}); continue
            if not out.exists():
                require(min(shutil.disk_usage(w).free,shutil.disk_usage(base).free)>=MIN_FREE,'Less than 8 GiB free; stopping safely. Free space and rerun.')
            try:
                status,r,source=process(j,climate,prov,workflow,utils,MonthlyTIModel)
                s=r['diagnostic_summary'];entry={'rgi_id':gid,'status':status,'rgi_area_km2':float(j['rgi_area_km2']),'modeled_area_km2':r['glacier_area_m2']/1e6,'melt_m3':s['total_disaggregated_melt_m3'],'rain_m3':s['total_disaggregated_rain_m3']}
                if a.prune_gdirs:
                    # Verify the published output before deleting the extracted source.
                    verify_existing(out,j,prov)
                    try: prune(source,j)
                    except Exception as e: entry['prune_warning']=str(e);print('PRUNE WARNING',gid,str(e),flush=True)
                consecutive=0
            except Exception as e:
                entry={'rgi_id':gid,'status':'BLOCKED','error_class':type(e).__name__,'error':str(e)};consecutive+=1
                with errors.open('a') as f: f.write(f'\n{utc()} {gid}\n{traceback.format_exc()}\n')
            entries.append(entry);print(f'[{i}/{EXPECTED}]',gid,entry['status'],entry.get('error',''),flush=True)
            if i%25==0: atomic_json(report,summary(entries,len(jobs),prov,False))
            require(consecutive<20,'20 consecutive failures; stopping to avoid repeated systemic errors. Inspect production_errors.log.')
        finished=True
    finally:
        atomic_json(report,summary(entries,len(jobs),prov,finished));print('Report:',report,flush=True)
    s=summary(entries,len(jobs),prov,finished);print(s['status'], 'complete',s['completed'],'blocked',s['blocked'],flush=True)
    if s['completed']!=EXPECTED: print('Incomplete inventory; export will require explicit accepted coverage.',flush=True)

def self_test():
    physics_self_test()
    ref = {'source': 'ERA5', 'monthly_raw_comparison': [{'month': m, 'temperature_difference_C': 0.001, 'precipitation_difference_mm': 0.08, 'oggm_uncalibrated_monthly_mm': 10.0} for m in MONTHS]}
    climate_gate(ref)
    ref['monthly_raw_comparison'][0]['temperature_difference_C'] = 1.0
    try:
        climate_gate(ref)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Climate mismatch accepted')
    s = summary([{'status': 'BLOCKED', 'rgi_id': 'test'}], EXPECTED, {}, False)
    assert s['completed'] == 0 and s['blocked'] == 1
    print('BASIN RUNNER SELF TEST PASS')

def physics_self_test():
    assert exact_index(np.array([38.0, 37.75]), 37.75) == 1
    try:
        exact_index(np.array([38.0, 37.75]), 37.8)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Off-grid coordinate accepted')
    fields = []
    for y in range(2005, 2020):
        for m in range(1, 13):
            n = calendar.monthrange(y, m)[1]
            one = np.ones((n, 1))
            fields.append({'days': [date(y, m, i + 1) for i in range(n)], 'pdd': one, 'snow': one, 'rain': one, 'monthly_melt_ref_per_bin': np.array([n * 2.0]), 'monthly_rain_ref_per_bin': np.array([float(n)]), 'monthly_snow_ref_per_bin': np.array([float(n)]), 'monthly_precip_ref_per_bin': np.array([n * 2.0]), 'precip_scaling': 1.0})
    daily, monthly, stats = physics.glacier_daily(fields, np.array([1000.0]))
    assert len(daily) == 5478 and len(monthly) == 180 and (stats['total_disaggregated_melt_m3'] == 10956)
    print('SELF TEST PASS: exact grid lookup and existing 180-month accounting kernel')


def provenance(native_hash, reproduce_omissions=False):
    return {'method':METHOD,'native_era5_sha256':native_hash,
            'worklist_sha256':sha(BASIN/'worklist.csv'),
            'physics_sha256':sha(Path(physics.__file__)),
            'driver_sha256':sha(Path(__file__)), 'archive_url':URL,
            'oggm_version':CONFIG['oggm_version'],
            'reproduce_omissions':bool(reproduce_omissions)}

if __name__=='__main__': main()
