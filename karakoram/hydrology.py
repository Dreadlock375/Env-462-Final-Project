"""Run CWatM, repeat calibration if requested, and evaluate frozen parameters."""
from __future__ import annotations
import argparse, csv, datetime as dt, fcntl, hashlib, json, math, os, re, subprocess, shutil
from pathlib import Path
import numpy as np
from .common import ROOT, CONFIG, BASIN, CLIMATE, MAPS, HYDROLOGY, OBS, SPATIAL, require, sha, atomic_json, read_json, read_csv, write_csv
FORCING=CLIMATE
MASK=BASIN/'mask.nc'
AREA=SPATIAL/'routing/cellarea.nc'
UNITS=dict(pr='kg m-2 s-1',tas='K',tasmin='K',tasmax='K',ps='Pa',huss='kg kg-1',sfcwind='m s-1',rsds='W m-2',rlds='W m-2')
def verify_forcing():
    import numpy as np
    from netCDF4 import Dataset, num2date
    report = json.loads((FORCING / 'conversion_report.json').read_text())
    require(report['status'] == 'DAILY_FORCING_PREPARED_NOT_RUNTIME_TESTED' and report['days'] == 5478, 'Unexpected conversion report')
    lat = np.arange(37.25, 30.74, -0.5)
    lon = np.arange(72.25, 82.26, 0.5)

    def crop(ds, name):
        iy = [int(np.argmin(abs(ds['lat'][:] - v))) for v in lat]
        ix = [int(np.argmin(abs(ds['lon'][:] - v))) for v in lon]
        require(np.array_equal(ds['lat'][:][iy], lat) and np.array_equal(ds['lon'][:][ix], lon), 'Static grid mismatch')
        return np.asarray(ds[name][:])[np.ix_(iy, ix)]
    with Dataset(MASK) as ds:
        mask = crop(ds, 'mask') == 1
    with Dataset(AREA) as ds:
        name = [v for v in ds.variables if ds[v].dimensions == ('lat', 'lon')]
        require(len(name) == 1, 'Ambiguous area variable')
        area = crop(ds, name[0]).astype('float64')
    require(mask.sum() == 67 and np.isfinite(area[mask]).all() and np.all(area[mask] > 0), 'Invalid basin mask or area')
    w = np.where(mask, area, 0)
    w /= w.sum()
    expected = [(d.year, d.month, d.day) for d in (dt.date(2005, 1, 1) + dt.timedelta(days=i) for i in range(5478))]
    verified = {}
    temps = {}
    basin = {}
    for var, units in UNITS.items():
        rel = f'cwatm_050deg/era5_{var}_daily_2005_2019.nc'
        path = FORCING / rel
        digest = sha(path)
        require(digest == report['output_sha256'][rel], f'Checksum mismatch: {path}')
        with Dataset(path) as ds:
            require(np.array_equal(ds['lat'][:], lat) and np.array_equal(ds['lon'][:], lon), f'Grid mismatch: {var}')
            tv = ds['time']
            dates = num2date(tv[:], tv.units, tv.calendar)
            require([(x.year, x.month, x.day) for x in dates] == expected and all((x.hour == x.minute == x.second == 0 for x in dates)), f'Date mismatch: {var}')
            v = ds[var]
            require(v.dimensions == ('time', 'lat', 'lon') and v.shape == (5478, 14, 21) and (v.units == units), f'Shape or units mismatch: {var}')
            ar = v[:]
            require(not np.ma.getmaskarray(ar).any() and np.isfinite(ar).all(), f'Invalid values: {var}')
            ar = np.asarray(ar, dtype='float64')
            if var in ('pr', 'rsds', 'rlds', 'sfcwind', 'huss'):
                require(ar.min() >= 0, f'Negative {var}')
            if var in ('tas', 'tasmin', 'tasmax'):
                temps[var] = ar
            if var in ('tas', 'pr'):
                basin[var] = np.sum(ar * w, axis=(1, 2))
            verified[var] = {'sha256': digest, 'min': float(ar.min()), 'max': float(ar.max()), 'units': units}
        print('Verified', var, flush=True)
    require(np.all(temps['tasmin'] <= temps['tas'] + 0.0001) and np.all(temps['tas'] <= temps['tasmax'] + 0.0001), 'Temperature extrema ordering error')
    annual = {}
    years = np.array([d[0] for d in expected])
    for year in range(2005, 2020):
        sel = years == year
        annual[str(year)] = {'precipitation_mm': float(basin['pr'][sel].sum() * 86400), 'temperature_mean_C': float(basin['tas'][sel].mean() - 273.15)}
        prior = report['annual_basin'][str(year)]
        require(abs(annual[str(year)]['precipitation_mm'] - prior['precipitation_mm']) < 0.01 and abs(annual[str(year)]['temperature_mean_C'] - prior['temperature_mean_C']) < 0.001, f'Output versus report discrepancy: {year}')
    return {'files': verified, 'annual_basin': annual, 'conversion_plan': report['plan']}

def array(values):
    a=np.asarray(values,dtype=float)
    require(np.isfinite(a).all() and (a>=0).all(),'Nonfinite or negative discharge')
    return a

def verify_reference_observations(values):
    """Verify the private reference observations without redistributing them."""
    reference = read_json(ROOT / 'reference/observation_fingerprint.json')
    encoded = np.asarray(array(values), dtype='<f8').tobytes()
    require(len(values) == reference['days'] and
            hashlib.sha256(encoded).hexdigest() == reference['sha256'],
            'Reference observations differ')


def metrics(s, o, dates):
    s, o = (array(s), array(o))
    require(len(s) == len(o) == len(dates) and len(s) > 1, 'Metric input length mismatch')
    r = float(np.corrcoef(s, o)[0, 1])
    alpha = float(s.std() / o.std())
    beta = float(s.mean() / o.mean())
    y = np.array([d.year for d in dates])
    m = np.array([d.month for d in dates])
    peaks = []
    for year in sorted(set(y)):
        ix = np.flatnonzero(y == year)
        i, j = (ix[np.argmax(s[ix])], ix[np.argmax(o[ix])])
        peaks.append(dict(year=int(year), observed_m3s=float(o[j]), simulated_m3s=float(s[i]), observed_date=str(dates[j]), simulated_date=str(dates[i]), relative_error_pct=float(100 * (s[i] / o[j] - 1)), timing_error_days=(dates[i] - dates[j]).days))
    sm = [float(s[m == i].mean()) for i in range(1, 13)]
    om = [float(o[m == i].mean()) for i in range(1, 13)]
    return dict(days=len(s), nse=float(1 - ((s - o) ** 2).sum() / ((o - o.mean()) ** 2).sum()), kge=float(1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)), r=r, std_ratio=alpha, bias_pct=100 * (beta - 1), rmse_m3s=float(np.sqrt(np.mean((s - o) ** 2))), mean_m3s=float(s.mean()), observed_mean_m3s=float(o.mean()), annual_peak_bias_pct=float(100 * (sum((p['simulated_m3s'] for p in peaks)) / sum((p['observed_m3s'] for p in peaks)) - 1)), annual_peak_abs_relative_error_pct=float(np.mean([abs(p['relative_error_pct']) for p in peaks])), annual_peak_mean_abs_timing_error_days=float(np.mean([abs(p['timing_error_days']) for p in peaks])), monthly_mean_abs_error_m3s=float(np.mean(np.abs(np.array(sm) - om))), monthly_means_m3s=sm, observed_monthly_means_m3s=om, annual_peaks=peaks)


CALIBRATION_END=dt.date.fromisoformat(CONFIG['calibration_end'])
FIRST=dt.date.fromisoformat(CONFIG['spinup_end'])+dt.timedelta(days=1)
LAST=dt.date.fromisoformat(CONFIG['end'])
FROZEN=CONFIG['frozen_parameters']


def dates_through(end):
    return [FIRST+dt.timedelta(days=i) for i in range((end-FIRST).days+1)]


def replace_setting(text,key,value):
    pattern=rf'(?m)^([ \t]*{re.escape(key)}[ \t]*=[ \t]*)([^\r\n]*)$'
    result,n=re.subn(pattern,lambda m:m[1]+str(value),text)
    require(n==1,f'Expected exactly one setting {key}, found {n}')
    return result


def settings(folder,end,params):
    text=(ROOT/'config/cwatm.ini.template').read_text().replace('@ROOT@',str(ROOT)).replace('@OUTPUT@',str(folder))
    text=replace_setting(text,'StepEnd',end.strftime('%d/%m/%Y'))
    for key,value in params.items():text=replace_setting(text,key,value)
    return text


def read_discharge(path, dates):
    records={}
    with path.open() as f:
        for row in csv.reader(f):
            if not row or not re.fullmatch(r'\d{2}/\d{2}/\d{4}',row[0].strip()):continue
            d=dt.datetime.strptime(row[0].strip(),'%d/%m/%Y').date()
            require(d not in records,f'Duplicate model date: {d}')
            records[d]=float(row[1])
    require(sorted(records)==dates,f'Model date coverage mismatch: {path}')
    return array([records[d] for d in dates])


def observations(dates):
    require(sha(OBS/'daily.csv')==read_json(OBS/'report.json')['daily_sha256'],'Prepared observations changed')
    records={dt.date.fromisoformat(r['date']):float(r['discharge_m3s']) for r in read_csv(OBS/'daily.csv')}
    return array([records[d] for d in dates])


def inputs():
    """Authenticate forcing files before running any hydrology case."""
    from .coupling import check_maps
    from .setup_inputs import static_hashes
    static=static_hashes(ROOT)
    source_manifest=ROOT/'config/source_manifest.json'
    if source_manifest.is_file():
        require(static==read_json(source_manifest)['static_sha256'],'Static map bytes differ from source manifest')
    local=ROOT/'config/local_inputs.json'
    if local.is_file():
        require(static==read_json(local)['static_sha256'],'Staged static map bytes changed')
    verified=verify_forcing()
    export=read_json(MAPS/'export_report.json')
    check_maps(MAPS,export)
    # The original model revision is part of the scientific configuration.
    repo=ROOT/'cwatm/source'
    revision=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
    require(revision==CONFIG['cwatm_revision'],f'CWatM revision differs: {revision}')
    status=subprocess.check_output(['git','-C',str(repo),'status','--porcelain','--untracked-files=no'],text=True).strip()
    require(not status,'CWatM tracked source modified')
    return dict(cwatm_revision=revision,static_sha256=static,climate_report=sha(CLIMATE/'conversion_report.json'),glacier_report=sha(MAPS/'export_report.json'),observation_report=sha(OBS/'report.json'),configuration=sha(ROOT/'config/project.json'),template=sha(ROOT/'config/cwatm.ini.template'))


def run_case(name, end, params, provenance):
    """Execute or reuse a provenance-matched model run; preserve interrupted outputs."""
    folder = HYDROLOGY / name
    folder.parent.mkdir(parents=True, exist_ok=True)
    text = settings(folder, end, params)
    days = dates_through(end)
    lock = folder.with_name(folder.name + '.lock')
    with lock.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if folder.exists():
            require((folder / 'settings.ini').is_file() and (folder / 'settings.ini').read_text() == text, 'Existing configuration differs')
            require(read_json(folder / 'provenance.json') == provenance, 'Existing provenance differs')
            if (folder / 'exit_code.txt').exists() and (folder / 'exit_code.txt').read_text().strip() == '0':
                series = read_discharge(folder / 'discharge_daily.csv', days)
                require(sha(folder / 'discharge_daily.csv') == read_json(folder / 'run_report.json')['discharge_sha256'], 'Existing discharge changed')
                print('Reused', name, flush=True)
                return series
            folder.rename(folder.with_name(folder.name + '.interrupted_' + dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')))
        folder.mkdir()
        (folder / 'settings.ini').write_text(text)
        atomic_json(folder / 'provenance.json', provenance)
        print('Running', name, 'log:', folder / 'run.log', flush=True)
        with (folder / 'run.log').open('w') as f:
            proc = subprocess.run(['bash', str(ROOT / 'run_cwatm.sh'), str(folder / 'settings.ini')], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env={**os.environ, 'PYTHONUNBUFFERED': '1'})
        require(proc.returncode == 0, f'CWatM failed: {folder}/run.log')
        series = read_discharge(folder / 'discharge_daily.csv', days)
        atomic_json(folder / 'run_report.json', dict(discharge_sha256=sha(folder / 'discharge_daily.csv'), settings_sha256=sha(folder / 'settings.ini'), days=len(series)))
        (folder / 'exit_code.txt').write_text('0\n')
        return series


SENSITIVITY={
    'precip_080':{'precipitation_coversion':'69.12'},
    'precip_070':{'precipitation_coversion':'60.48'},
    'crop_125':{'crop_correct':'1.25'},
    'snow_075':{'SnowMeltCoef':'0.003'},
    'summer_snow_050':{'IceMeltCoef':'0.0035'},
    'routing_200':{'manningsN':'2.0'},
    'concentration_200':{'runoffConc_factor':'2.0'},
}
P70=SENSITIVITY['precip_070']
COMBINED={
    'p70_s50':{**P70,'IceMeltCoef':'0.0035'},
    'p70_r15':{**P70,'manningsN':'1.5'},
    'p70_r20':{**P70,'manningsN':'2.0'},
    'p70_s50_r15':{**P70,'IceMeltCoef':'0.0035','manningsN':'1.5'},
    'p70_s50_r20':{**P70,'IceMeltCoef':'0.0035','manningsN':'2.0'},
    'p70_s75_r15':{**P70,'IceMeltCoef':'0.00525','manningsN':'1.5'},
}


def calibrate(provenance):
    """Repeat the original seven sensitivity and six combined cases; use only calibration dates."""
    days = dates_through(CALIBRATION_END)
    o = observations(days)
    results = {}
    series = {}
    for name, params in {**SENSITIVITY, **COMBINED}.items():
        s = run_case('calibration/' + name, CALIBRATION_END, params, provenance)
        series[name] = s
        results[name] = metrics(s, o, days)
        print(name, results[name]['nse'], results[name]['annual_peak_abs_relative_error_pct'], flush=True)
    candidates = ['precip_070', *COMBINED]
    selected = min(candidates, key=lambda k: results[k]['annual_peak_abs_relative_error_pct'])
    report = dict(selection=selected, parameters={**SENSITIVITY, **COMBINED}[selected], cases=results, selection_basis='Minimum mean absolute relative annual peak magnitude error among precipitation-adjusted candidates, calibration period only')
    folder = HYDROLOGY / 'calibration'
    atomic_json(folder / 'report.json', report)
    write_csv(folder / 'daily.csv', [dict(date=str(d), observed_m3s=float(o[i]), **{k: float(v[i]) for k, v in series.items()}) for i, d in enumerate(days)])
    require(report['parameters'] == FROZEN, 'Calibration selected different parameters; inspect results instead of silently changing the historical reference')
    return series[selected]


def evaluate(provenance, repeat_calibration=False):
    """Replay baseline and frozen validation without any old diagnostic output."""
    days = dates_through(LAST)
    o = observations(days)
    baseline = run_case('uncalibrated', LAST, {}, provenance)
    cal = calibrate(provenance) if repeat_calibration else None
    s = run_case('calibrated', LAST, FROZEN, provenance)
    cut = len(dates_through(CALIBRATION_END))
    if cal is not None:
        require(np.allclose(s[:cut], cal, rtol=1e-08, atol=1e-05), 'Calibration-period reproduction failed')
    folder = HYDROLOGY / 'calibrated'
    periods = {}
    for name, sl in [('calibration_2006_2012', slice(0, cut)), ('validation_2013_2019', slice(cut, None)), ('full_2006_2019', slice(None))]:
        periods[name] = dict(uncalibrated=metrics(baseline[sl], o[sl], days[sl]), calibrated=metrics(s[sl], o[sl], days[sl]))
    result = dict(status='FROZEN_PARAMETER_EVALUATION_COMPLETE', parameters=FROZEN, periods=periods, provenance=provenance, calibration_search_repeated=repeat_calibration)
    atomic_json(folder / 'validation_report.json', result)
    write_csv(folder / 'validation_daily.csv', [dict(date=str(d), period='calibration' if i < cut else 'validation', observed_m3s=float(o[i]), uncalibrated_m3s=float(baseline[i]), calibrated_m3s=float(s[i])) for i, d in enumerate(days)])
    reference = read_csv(ROOT / 'reference/validation_daily.csv')
    require([r['date'] for r in reference] == [str(d) for d in days], 'Reference dates mismatch')
    verify_reference_observations(o)
    orig = array([r['calibrated_m3s'] for r in reference])
    diff = s - orig
    regression = dict(daily_match=bool(np.allclose(s, orig, rtol=1e-08, atol=1e-05)), max_abs_m3s=float(np.max(np.abs(diff))), rmse_difference_m3s=float(np.sqrt(np.mean(diff ** 2))))
    atomic_json(folder / 'regression_report.json', regression)
    print('Regression:', regression, flush=True)
    for name, values in periods.items():
        print(name, {k: values['calibrated'][k] for k in ('nse', 'kge', 'bias_pct', 'annual_peak_abs_relative_error_pct')}, flush=True)
    return regression


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeat-calibration',action='store_true')
    parser.add_argument('--require-reference-match',action='store_true',help='Fail if daily calibrated series does not match original reference')
    args=parser.parse_args()
    result=evaluate(inputs(),args.repeat_calibration)
    if args.require_reference_match:require(result['daily_match'],'Full historical daily regression failed')

if __name__=='__main__':main()
