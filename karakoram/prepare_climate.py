#!/usr/bin/env python3
"""Merge authenticated CDS GRIB and ECMWF ARCO tiles, produce UTC daily forcing.
Run from project root. No network access, no model settings or raw files changed.
"""
import argparse, calendar, collections, datetime as dt, hashlib, json, re, shutil, tempfile
from pathlib import Path
import numpy as np

from .common import ROOT, BASIN
RAW=ROOT/'data/climate/era5/raw'
OUT=ROOT/'data/climate/era5/processed/daily_2005_2019'
MASK=BASIN/'mask.nc'
AREA=ROOT/'data/cwatm/spatial/30min/cwatm_input_30min/routing/cellarea.nc'
START=np.datetime64('2005-01-01T00','h'); END=np.datetime64('2020-01-01T00','h')
TIMES=np.arange(START,END+np.timedelta64(1,'h'),np.timedelta64(1,'h'))
LAT=np.arange(38.,29.99,-.25); LON=np.arange(71.5,83.01,.25)
TLAT=np.arange(37.25,30.74,-.5); TLON=np.arange(72.25,82.26,.5)
VARS=('t2m','d2m','sp','u10','v10','tp','ssrd','strd')
ACC={'tp','ssrd','strd'}
# Explicit tolerances in the original accumulation units.
NEG_TOL={'tp':1e-7,'ssrd':10.,'strd':1.}
CLIPPED={v:{'count':0,'minimum':0.,'sum_removed':0.,'maximum_daily_cell_correction':0.} for v in ACC}
ALIAS={'2t':'t2m','2d':'d2m','10u':'u10','10v':'v10'}
META={'pr':('kg m-2 s-1','precipitation_flux'), 'tas':('K','air_temperature'), 'tasmin':('K','air_temperature'), 'tasmax':('K','air_temperature'), 'ps':('Pa','surface_air_pressure'), 'huss':('kg kg-1','specific_humidity'), 'sfcwind':('m s-1','wind_speed'), 'rsds':('W m-2','surface_downwelling_shortwave_flux_in_air'), 'rlds':('W m-2','surface_downwelling_longwave_flux_in_air')}

def require(ok,msg):
    if not ok: raise ValueError(msg)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def exact_indices(full, sub):
    full = np.asarray(full)
    sub = np.asarray(sub)
    indices = []
    for x in sub:
        q = np.flatnonzero(np.isclose(full, x, rtol=0, atol=1e-07))
        require(len(q) == 1, f'Coordinate {x} not unique in expected grid')
        indices.append(q[0])
    require(len(set(indices)) == len(indices), 'Duplicate coordinates')
    return np.asarray(indices, dtype=int)

def weights(src,dst,lat=False):
    # Rectangular cell overlap, latitude weighted by spherical area.
    lo=np.maximum(dst[:,None]-.25,src[None,:]-.125)
    hi=np.minimum(dst[:,None]+.25,src[None,:]+.125)
    if lat:
        w=np.where(hi>lo,np.sin(np.deg2rad(hi))-np.sin(np.deg2rad(lo)),0)
        den=np.sin(np.deg2rad(dst+.25))-np.sin(np.deg2rad(dst-.25))
    else:
        w=np.maximum(hi-lo,0); den=np.full(len(dst),.5)
    w=w/den[:,None]
    require(np.allclose(w.sum(axis=1),1,rtol=0,atol=1e-12),'Incomplete conservative remap coverage')
    return w

WY=weights(LAT,TLAT,True); WX=weights(LON,TLON)

def remap(x): return WY @ x @ WX.T

def daily(a):
    # Instantaneous hours 00..23; preceding-hour accumulations 01..next 00.
    t=a['t2m'][:24].astype('float64'); d=a['d2m'][:24].astype('float64'); p=a['sp'][:24].astype('float64')
    require(np.all((t>170)&(t<350)&(d>160)&(d<350)&(p>25000)&(p<110000)), 'Implausible hourly temperature or pressure')
    # Saturation vapour pressure over water (ECMWF dew point convention), Pa.
    e=611.21*np.exp(17.502*(d-273.16)/(d-32.19))
    require(np.all(e<p),'Vapour pressure exceeds surface pressure')
    q=.622*e/(p-.378*e)
    wind=np.hypot(a['u10'][:24].astype('float64'),a['v10'][:24].astype('float64'))
    sums={}
    for v in ACC:
        x=a[v][1:25].astype('float64')
        # Reject material negative accumulations. Only tiny rounding noise may be clipped.
        tol=NEG_TOL[v]
        minimum=float(x.min())
        require(minimum>=-tol,f'Material negative accumulation: {v}: {minimum} (tolerance {tol})')
        negative=x<0
        if negative.any():
            stat=CLIPPED[v]; removed=np.maximum(-x,0)
            stat['count']+=int(negative.sum())
            stat['minimum']=min(stat['minimum'],minimum)
            stat['sum_removed']+=float(removed.sum())
            stat['maximum_daily_cell_correction']=max(stat['maximum_daily_cell_correction'],float(removed.sum(axis=0).max()))
        sums[v]=np.maximum(x,0).sum(axis=0)
    return {'tas':t.mean(0),'tasmin':t.min(0),'tasmax':t.max(0),'ps':p.mean(0),'huss':q.mean(0),'sfcwind':wind.mean(0),'pr':sums['tp']*1000/86400,'rsds':sums['ssrd']/86400,'rlds':sums['strd']/86400}

def month_span(y,m):
    start=np.datetime64(f'{y:04d}-{m:02d}','M')
    return start.astype('datetime64[h]'),(start+1).astype('datetime64[h]')

def manifest_plan(path):
    doc = json.loads(path.read_text())
    ident = doc['plan']
    require(path.parent.name == ident, 'Manifest directory does not match plan ID')
    require(set(doc['variables']) == set(VARS) and doc['spatial_points'] == len(LAT) * len(LON), 'Unexpected manifest grid/variables')
    covered = np.zeros(len(TIMES), bool)
    cds = []
    for item in doc['reused_cds']:
        a, b = month_span(item['year'], item['month'])
        mask = TIMES == a if item['boundary'] else (TIMES >= a) & (TIMES < b)
        require(mask.any() and (not np.any(covered & mask)), 'CDS overlap/out of period')
        covered |= mask
        p = ROOT / item['file']
        require(p.is_file() and sha(p) == item['sha256'], f'CDS checksum failed: {p}')
        cds.append((p, TIMES[mask]))
    missing = TIMES[~covered]
    require(len(missing) == doc['hours'], 'ARCO hour count mismatch')
    cover = {v: np.zeros((len(LAT), len(LON)), np.int16) for v in VARS}
    tiles = []
    from netCDF4 import Dataset, date2num
    for k, item in enumerate(doc['tiles'], 1):
        v = item['variable']
        p = path.parent / item['file']
        require(item['plan'] == ident and item['hours'] == len(missing), f'Tile plan mismatch: {p}')
        require(p.is_file() and p.stat().st_size == item['bytes'] and (sha(p) == item['sha256']), f'Tile checksum failed: {p}')
        with Dataset(p) as ds:
            f = ds[v]
            require(f.dimensions == ('time', 'latitude', 'longitude'), f'Unexpected tile dimensions {p}')
            iy = exact_indices(LAT, ds['latitude'][:])
            ix = exact_indices(LON, ds['longitude'][:])
            require(f.shape == (len(missing), len(iy), len(ix)), f'Tile shape mismatch {p}')
            tv = ds['time']
            units = tv.units
            cal = getattr(tv, 'calendar', 'standard')
            base = dt.datetime(2005, 1, 1)
            offset = date2num(base, units, cal)
            scale = date2num(base + dt.timedelta(hours=1), units, cal) - offset
            expected = offset + (missing - START).astype('int64') * scale
            require(np.allclose(tv[:], expected, rtol=0, atol=abs(scale) * 1e-05), f'Tile time mismatch {p}')
            require(getattr(ds, 'acquisition_plan', '') == ident, f'Tile internal plan mismatch {p}')
            cover[v][np.ix_(iy, ix)] += 1
        tiles.append((v, p, iy, ix))
        if k % 100 == 0 or k == len(doc['tiles']):
            print(f"Verified {k}/{len(doc['tiles'])} tiles", flush=True)
    if len(missing):
        require(all((np.all(x == 1) for x in cover.values())), 'Tile spatial gap or overlap')
    print(f'Authenticated {len(cds)} CDS files and {len(tiles)} ARCO tiles; {len(TIMES)} planned hours', flush=True)
    return (doc, cds, missing, tiles)

def read_grib(path, expected, wanted, arrays, start):
    from eccodes import codes_grib_new_from_file, codes_get, codes_get_values, codes_release
    seen = {v: set() for v in VARS}
    valid = set(expected.astype('int64'))
    selected = set(wanted.astype('int64'))
    with path.open('rb') as f:
        while True:
            g = codes_grib_new_from_file(f)
            if g is None:
                break
            try:
                raw = codes_get(g, 'shortName')
                v = ALIAS.get(raw, raw)
                require(v in seen, f'Unexpected GRIB variable {raw}')
                date = str(codes_get(g, 'validityDate'))
                hm = int(codes_get(g, 'validityTime'))
                require(hm % 100 == 0, f'Non-hourly GRIB time: {path}')
                h = np.datetime64(f'{date[:4]}-{date[4:6]}-{date[6:]}T{hm // 100:02d}', 'h')
                stamp = int(h.astype('int64'))
                require(stamp in valid and stamp not in seen[v], f'Unexpected/duplicate GRIB hour {path}: {v} {h}')
                seen[v].add(stamp)
                require((codes_get(g, 'Ni'), codes_get(g, 'Nj')) == (47, 33), f'GRIB grid shape {path}')
                for key, val in [('latitudeOfFirstGridPointInDegrees', 38), ('latitudeOfLastGridPointInDegrees', 30), ('longitudeOfFirstGridPointInDegrees', 71.5), ('longitudeOfLastGridPointInDegrees', 83), ('iDirectionIncrementInDegrees', 0.25), ('jDirectionIncrementInDegrees', 0.25)]:
                    require(abs(float(codes_get(g, key)) - val) < 1e-06, f'GRIB grid {key}: {path}')
                require(int(codes_get(g, 'scanningMode')) == 0, f'Unexpected GRIB scan order: {path}')
                require(codes_get(g, 'stepType') == ('accum' if v in ACC else 'instant'), f'Unexpected step type: {path} {v}')
                if v in ACC:
                    require(str(codes_get(g, 'stepUnits')) in ('1', 'h') and float(codes_get(g, 'endStep')) - float(codes_get(g, 'startStep')) == 1, f'Not preceding one-hour accumulation: {path} {v}')
                if stamp in selected:
                    a = np.asarray(codes_get_values(g)).reshape(33, 47)
                    require(int(codes_get(g, 'numberOfMissing')) == 0 and np.isfinite(a).all(), f'Missing GRIB data {path}')
                    arrays[v][int((h - start).astype('int64'))] = a
            finally:
                codes_release(g)
    require(all((s == valid for s in seen.values())), f'Incomplete GRIB record coverage: {path}')

def load_year(y, cds, missing, tiles):
    from netCDF4 import Dataset
    a = np.datetime64(f'{y}-01-01T00', 'h')
    b = np.datetime64(f'{y + 1}-01-01T00', 'h')
    count = int((b - a).astype('int64')) + 1
    arrays = {v: np.full((count, 33, 47), np.nan, dtype='float32') for v in VARS}
    for p, expected in cds:
        wanted = expected[(expected >= a) & (expected <= b)]
        if len(wanted):
            print('  CDS', p.name, flush=True)
            read_grib(p, expected, wanted, arrays, a)
    lo = int(np.searchsorted(missing, a))
    hi = int(np.searchsorted(missing, b, side='right'))
    if hi > lo:
        indices = (missing[lo:hi] - a).astype('int64')
        for j, (v, p, iy, ix) in enumerate(tiles, 1):
            with Dataset(p) as ds:
                raw = ds[v][lo:hi, :, :]
                require(not np.ma.getmaskarray(raw).any(), f'Masked tile values: {p}')
                x = np.asarray(raw, dtype='float32')
                require(np.isfinite(x).all(), f'Non-finite tile data: {p}')
                arrays[v][np.ix_(indices, iy, ix)] = x
            if j % 234 == 0 or j == len(tiles):
                print(f'  ARCO {j}/{len(tiles)} tiles', flush=True)
    for v, x in arrays.items():
        require(np.isfinite(x).all(), f'Missing/unfilled hours for {y} {v}')
    return arrays

def create_nc(path, lat, lon, variables):
    from netCDF4 import Dataset
    f = Dataset(path, 'w', format='NETCDF4')
    f.createDimension('time', None)
    f.createDimension('lat', len(lat))
    f.createDimension('lon', len(lon))
    t = f.createVariable('time', 'f8', ('time',))
    t.units = 'days since 1900-01-01 00:00:00'
    t.calendar = 'proleptic_gregorian'
    t.standard_name = 'time'
    t.axis = 'T'
    for k, data, unit, axis in [('lat', lat, 'degrees_north', 'Y'), ('lon', lon, 'degrees_east', 'X')]:
        x = f.createVariable(k, 'f8', (k,))
        x[:] = data
        x.units = unit
        x.standard_name = {'lat': 'latitude', 'lon': 'longitude'}[k]
        x.axis = axis
    for v in variables:
        x = f.createVariable(v, 'f4', ('time', 'lat', 'lon'), zlib=True, complevel=2, chunksizes=(min(31, 5478), len(lat), len(lon)), fill_value=False)
        x.units, x.standard_name = META[v]
        x.cell_methods = 'time: ' + ('minimum' if v == 'tasmin' else 'maximum' if v == 'tasmax' else 'mean')
    f.source = 'ECMWF ERA5; CDS GRIB and ECMWF ARCO Zarr'
    f.history = 'UTC daily aggregation; preceding-hour accumulations 01 UTC to next day 00 UTC'
    f.spatial_method = 'native 0.25 degree' if len(lat) == 33 else 'spherical rectangular area-overlap weighted mean of daily native fields; no bias correction'
    return f

def basin_weights():
    from netCDF4 import Dataset
    with Dataset(MASK) as f:
        iy = exact_indices(f['lat'][:], TLAT)
        ix = exact_indices(f['lon'][:], TLON)
        mask = np.asarray(f['mask'][:])[np.ix_(iy, ix)] == 1
    require(mask.sum() == 67, 'Unexpected basin mask cell count')
    with Dataset(AREA) as f:
        iy = exact_indices(f['lat'][:], TLAT)
        ix = exact_indices(f['lon'][:], TLON)
        var = [v for v in f.variables if f[v].dimensions == ('lat', 'lon')]
        require(len(var) == 1, 'Ambiguous cell area variable')
        area = np.asarray(f[var[0]][:])[np.ix_(iy, ix)]
    require(np.isfinite(area[mask]).all() and np.all(area[mask] > 0), 'Invalid basin cell area')
    w = np.where(mask, area, 0)
    return (w / w.sum(), float(w.sum() / 1000000.0))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, default=None)
    p.add_argument('--self-test', action='store_true')
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    if OUT.exists():
        report = json.loads((OUT / 'conversion_report.json').read_text())
        require(all((sha(OUT / n) == h for n, h in report['output_sha256'].items())), 'Existing climate output changed')
        print('Verified existing daily ERA5 forcing:', OUT)
        return
    if args.manifest is None:
        choices = sorted((RAW / 'ecmwf_arco').glob('*/COMPLETE.json'))
        require(len(choices) == 1, 'Specify --manifest: zero or multiple completed acquisition plans')
        args.manifest = choices[0]
    from netCDF4 import Dataset, date2num
    import eccodes
    doc, cds, missing, tiles = manifest_plan(args.manifest)
    bw, area = basin_weights()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.daily_era5_', dir=OUT.parent))
    print('Staging:', stage, flush=True)
    native = stage / 'native_025deg'
    native.mkdir()
    coarse = stage / 'cwatm_050deg'
    coarse.mkdir()
    handles = []
    monthly = {}
    ranges = {v: [float('inf'), float('-inf')] for v in META}
    count = 0
    try:
        nf = create_nc(native / 'era5_daily_2005_2019.nc', LAT, LON, META)
        handles.append(nf)
        cf = {}
        for v in META:
            cf[v] = create_nc(coarse / f'era5_{v}_daily_2005_2019.nc', TLAT, TLON, [v])
            handles.append(cf[v])
        for y in range(2005, 2020):
            print(f'Processing {y}', flush=True)
            a = load_year(y, cds, missing, tiles)
            nd = calendar.isleap(y) + 365
            for day in range(nd):
                date = dt.datetime(y, 1, 1) + dt.timedelta(days=day)
                idx = day * 24
                try:
                    d = daily({v: x[idx:idx + 25] for v, x in a.items()})
                except ValueError as exc:
                    raise ValueError(f'{date:%Y-%m-%d}: {exc}') from exc
                out = {v: remap(x) for v, x in d.items()}
                require(np.all(out['tasmin'] <= out['tas'] + 1e-07) and np.all(out['tas'] <= out['tasmax'] + 1e-07), 'Daily temperature ordering')
                tm = date2num(date, nf['time'].units, nf['time'].calendar)
                nf['time'][count] = tm
                for v in META:
                    nf[v][count] = d[v]
                    cf[v]['time'][count] = tm
                    cf[v][v][count] = out[v]
                    ranges[v][0] = min(ranges[v][0], float(out[v].min()))
                    ranges[v][1] = max(ranges[v][1], float(out[v].max()))
                m = monthly.setdefault(date.strftime('%Y-%m'), {'days': 0, 'precipitation_mm': 0.0, 'temperature_sum_C': 0.0})
                m['days'] += 1
                m['precipitation_mm'] += float(np.sum(out['pr'] * bw) * 86400)
                m['temperature_sum_C'] += float(np.sum((out['tas'] - 273.15) * bw))
                count += 1
            del a
            for f in handles:
                f.sync()
            print(f'  Daily records written: {count}', flush=True)
        require(count == 5478, 'Wrong daily record count')
    finally:
        for f in handles:
            f.close()
    for f in [native / 'era5_daily_2005_2019.nc', *sorted(coarse.glob('*.nc'))]:
        with Dataset(f) as ds:
            require(len(ds.dimensions['time']) == 5478 and np.all(np.diff(ds['time'][:]) == 1), f'Output time axis failed {f}')
            require(float(ds['time'][0]) == date2num(dt.datetime(2005, 1, 1), ds['time'].units, ds['time'].calendar) and float(ds['time'][-1]) == date2num(dt.datetime(2019, 12, 31), ds['time'].units, ds['time'].calendar), f'Output dates failed {f}')
            for v in META.keys() & ds.variables.keys():
                for i in range(0, 5478, 366):
                    x = ds[v][i:i + 366]
                    require(not np.ma.getmaskarray(x).any() and np.isfinite(x).all(), f'Output invalid {f} {v}')
    for m in monthly.values():
        m['temperature_mean_C'] = m.pop('temperature_sum_C') / m['days']
    annual = {}
    for key, m in monthly.items():
        o = annual.setdefault(key[:4], {'days': 0, 'precipitation_mm': 0.0, 'temperature_sum_C': 0.0})
        o['days'] += m['days']
        o['precipitation_mm'] += m['precipitation_mm']
        o['temperature_sum_C'] += m['temperature_mean_C'] * m['days']
    for o in annual.values():
        o['temperature_mean_C'] = o.pop('temperature_sum_C') / o['days']
    report = {'status': 'DAILY_FORCING_PREPARED_NOT_RUNTIME_TESTED', 'plan': doc['plan'], 'days': count, 'cds_files': len(cds), 'arco_tiles': len(tiles), 'native_grid': [33, 47], 'cwatm_grid': [14, 21], 'basin_cells': 67, 'basin_area_km2': area, 'daily_ranges_cwatm': ranges, 'monthly_basin': monthly, 'annual_basin': annual, 'output_sha256': {str(f.relative_to(stage)): sha(f) for f in stage.rglob('*.nc')}, 'method': {'utc_instant': '00-23 UTC', 'utc_accum': '01 UTC through following 00 UTC', 'precipitation': 'm accumulation to kg m-2 s-1', 'radiation': 'J m-2 accumulation /86400', 'wind': 'mean hourly vector magnitude at 10m', 'humidity': 'specific humidity from hourly dew point and pressure then daily mean', 'remap': 'spherical overlap; daily local extrema averaged spatially; no orographic correction'}}
    report['negative_accumulation_qc'] = {v: {**CLIPPED[v], 'absolute_tolerance': NEG_TOL[v], 'units': 'm' if v == 'tp' else 'J m-2'} for v in sorted(ACC)}
    print('Negative accumulation QC:', json.dumps(report['negative_accumulation_qc']), flush=True)
    (stage / 'conversion_report.json').write_text(json.dumps(report, indent=2) + '\n')
    stage.rename(OUT)
    print('DAILY FORCING PREPARED:', OUT, flush=True)
    print('Report:', OUT / 'conversion_report.json', flush=True)

def self_test():
    a = {v: np.zeros((25, 33, 47)) for v in VARS}
    a['t2m'][:] = 280
    a['d2m'][:] = 275
    a['sp'][:] = 80000
    a['u10'][:] = 3
    a['v10'][:] = 4
    a['tp'][:] = 0.001
    a['tp'][0] = 99
    a['tp'][24] = 0.002
    a['ssrd'][:] = 3600
    a['strd'][:] = 7200
    o = daily(a)
    require(np.allclose(o['pr'] * 86400, 25), 'Accumulation boundary self-test')
    require(np.allclose(o['rsds'], 1) and np.allclose(o['rlds'], 2), 'Radiation test')
    require(np.allclose(o['sfcwind'], 5), 'Wind test')
    require(np.allclose(remap(np.ones((33, 47))), 1), 'Remap constant test')
    require(np.allclose(remap(o['pr']) * 86400, 25), 'Remap precipitation test')
    print('Self-test PASS')
    a['ssrd'][1, 0, 0] = -1.1655092239379883
    o = daily(a)
    require(np.isclose(o['rsds'][0, 0], 23 / 24), 'Tiny negative radiation clipping test')
    require(CLIPPED['ssrd']['count'] == 1 and np.isclose(CLIPPED['ssrd']['sum_removed'], 1.1655092239379883), 'Clipping accounting test')
    a['ssrd'][1, 0, 0] = -11.0
    try:
        daily(a)
    except ValueError as exc:
        require('Material negative accumulation: ssrd' in str(exc), 'Unexpected negative test error')
    else:
        raise AssertionError('Large negative radiation was accepted')
    print('Negative accumulation QC tests PASS')
if __name__=='__main__': main()
