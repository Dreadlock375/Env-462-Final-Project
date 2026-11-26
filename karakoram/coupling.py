#!/usr/bin/env python3
"""Export ERA5 glacier fluxes and run a matched provisional historical CWatM experiment."""
from __future__ import annotations
import argparse, csv, fcntl, json, math, os, shutil, subprocess, tempfile, zipfile
from pathlib import Path
import numpy as np
from . import glaciers as p
from . import basin as spatial
from .common import ROOT, BASIN, MAPS, CONFIG, atomic_json, read_csv as read_table

PROD = ROOT / p.OUT
LAT = np.arange(37.25,30.74,-.5)
LON = np.arange(72.25,82.26,.5)
N = len(p.DAYS)
DATES = np.array([d.isoformat() for d in p.DAYS], dtype='U10')
NAMES = {'melt':'melt_on_glacier_2005_2019.nc','rain':'liq_prcp_on_glacier_2005_2019.nc','cover':'glacier_cover_fraction_2005_2019.nc'}
VARIABLES = {'melt':'melt_on_glacier','rain':'liq_prcp_on_glacier','cover':'on_area'}
require = p.require

def readj(path): return json.loads(path.read_text())
def writej(path, obj): atomic_json(path,obj)
def close(a,b): return np.allclose(a,b,rtol=1e-10,atol=0.01)

def regional_index(r,c):
    y,x = r-105,c-504
    require(1<=y<13 and 1<=x<20,f'Cell outside audited CWatM cut: {r},{c}')
    return y,x

def add_glacier(melt, rain, cover, job, overlap, route, values):
    gid = job['rgi_id']
    fraction = min(1.0, float(job['polygon_overlap_besham_km2']) / float(job['rgi_area_km2']))
    require(0 < fraction <= 1, f'Invalid overlap fraction {gid}')
    y, x = regional_index(*route)
    melt[:, y, x] += values[:, 0] * fraction
    rain[:, y, x] += values[:, 1] * fraction
    for area, r, c in overlap:
        y, x = regional_index(r, c)
        cover[y, x] += area
    return fraction

def read_daily(path):
    # Validate dates and numeric components without assuming column order.
    with path.open(newline='') as f:
        header=next(csv.reader(f)); require(header==['date','snowfall_m3','rain_m3','snowmelt_m3','firnmelt_m3','icemelt_m3','total_melt_m3','surface_liquid_m3','snow_firn_storage_m3_end','firn_to_ice_m3'],f'Unexpected header {path}')
        arr=np.loadtxt(f,delimiter=',',dtype=[('date','U10')]+[(x,'f8') for x in header[1:]],ndmin=1)
    require(arr.shape==(N,) and np.array_equal(arr['date'],DATES),f'Date mismatch {path}')
    for k in header[1:]:require(np.isfinite(arr[k]).all() and (arr[k]>=-1e-8).all(),f'Invalid {k}: {path}')
    require(close(arr['total_melt_m3'],arr['snowmelt_m3']+arr['firnmelt_m3']+arr['icemelt_m3']) and close(arr['surface_liquid_m3'],arr['total_melt_m3']+arr['rain_m3']),f'Component mismatch {path}')
    return np.maximum(np.column_stack((arr['total_melt_m3'],arr['rain_m3'])),0.)

def write_regional(path, var, data, meteo):
    from netCDF4 import Dataset
    with Dataset(meteo) as src, Dataset(path, 'w', format='NETCDF4') as dst:
        for name, n in [('time', N), ('lat', 14), ('lon', 21)]:
            dst.createDimension(name, n)
            sv = src[name]
            v = dst.createVariable(name, 'f8', (name,))
            v.setncatts({k: sv.getncattr(k) for k in sv.ncattrs() if k != '_FillValue'})
            v[:] = sv[:]
        v = dst.createVariable(var, 'f4', ('time', 'lat', 'lon'), zlib=True, complevel=3, shuffle=True, chunksizes=(31, 14, 21))
        v.units = 'm3 day-1'
        v[:] = data
        dst.title = 'ERA5 glacier potential liquid water; provisional overlap-scaled release mapping'
        dst.validation_status = 'NOT INDEPENDENTLY VALIDATED'

def check_maps(folder, report):
    from netCDF4 import Dataset
    for k, name in NAMES.items():
        path = folder / name
        require(p.sha(path) == report['output_sha256'][name], f'Map checksum mismatch {name}')
        with Dataset(path) as d:
            v = d[VARIABLES[k]]
            require(v.dimensions == ('time', 'lat', 'lon'), f'Dimensions {name}')
            if k == 'cover':
                require(v.shape == (15, 360, 720) and v.units == '1', 'Cover shape/units')
                require(np.array_equal(d['lat'][:], 89.75 - np.arange(360) * 0.5) and np.array_equal(d['lon'][:], -179.75 + np.arange(720) * 0.5), 'Cover grid')
                a = p.finite(v[:])
                require(a.min() >= 0 and a.max() <= 1 and np.all(a == a[0]), 'Invalid static cover')
                require(p.dates(d['time']) == [p.date(y, 1, 1) for y in range(2005, 2020)], 'Cover dates')
                require(np.count_nonzero(a[0]) == np.count_nonzero(a[0, 106:118, 505:524]), 'Cover outside model cut')
            else:
                require(v.shape == (N, 14, 21) and v.units == 'm3 day-1', 'Flux shape/units')
                require(np.array_equal(d['lat'][:], LAT) and np.array_equal(d['lon'][:], LON), 'Flux grid')
                require(p.dates(d['time']) == p.DAYS, 'Flux dates')
                a = p.finite(v[:])
                require(a.min() >= 0 and np.count_nonzero(a) == np.count_nonzero(a[:, 1:13, 1:20]), 'Flux outside model cut or negative')
                require(np.isclose(a.sum(), report['exported_m3'][k], rtol=1e-07, atol=1), 'Flux readback volume')

def export(accept_omissions=False):
    from netCDF4 import Dataset
    jobs = {j['rgi_id']: j for j in spatial.read_jobs(ROOT)}
    overlap = load_overlap(ROOT)
    cells, _, _, _ = load_grid(ROOT)
    require(set(jobs) == set(overlap), 'Worklist / overlap IDs differ')
    report_path = PROD / 'production_report.json'
    report = readj(report_path)
    rh = p.sha(report_path)
    require(report['expected_glaciers'] == report['selected_glaciers'] == len(jobs) == 13952 and report['not_processed_in_selection'] == 0, 'Production not fully attempted')
    entries = {e['rgi_id']: e for e in report['entries']}
    require(len(entries) == len(report['entries']) == len(jobs) and set(entries) == set(jobs), 'Production entry mismatch')
    ok = {gid for gid, e in entries.items() if e['status'] in ('COMPUTED', 'REUSED')}
    bad = set(jobs) - ok
    require(len(ok) == report['completed'] and len(bad) == report['blocked'] <= 7 and (len(ok) >= 13945) and (bad == {e['rgi_id'] for e in report['blocked_entries']}) and all((entries[g]['status'] == 'BLOCKED' for g in bad)), 'Expected at least 13,945 completed, at most seven blocked')
    area_total = math.fsum((float(j['polygon_overlap_besham_km2']) for j in jobs.values()))
    missing_area = math.fsum((float(jobs[g]['polygon_overlap_besham_km2']) for g in bad))
    require(missing_area <= CONFIG['max_missing_overlap_km2'], 'Omitted area exceeds accepted coverage limit')
    require(not bad or accept_omissions, 'Incomplete production: explicitly use --accept-omissions after reviewing production_report.json')
    meteo = ROOT / p.CLIMATE / 'cwatm_050deg/era5_pr_daily_2005_2019.nc'
    cr = readj(ROOT / p.CLIMATE / 'conversion_report.json')
    prov = report['provenance']
    expected = p.provenance(cr['output_sha256'][str(p.NATIVE)], prov['reproduce_omissions'])
    require(prov == expected, 'Production provenance mismatch')
    require(p.sha(meteo) == cr['output_sha256']['cwatm_050deg/' + meteo.name], 'ERA5 precipitation checksum')
    with Dataset(meteo) as ds:
        require(np.array_equal(ds['lat'][:], LAT) and np.array_equal(ds['lon'][:], LON) and (p.dates(ds['time']) == p.DAYS), 'ERA5 grid or dates')
    input_hashes = {str(x): p.sha(ROOT / x) for x in (OVERLAP, CELLS, ROUTES, WORKLIST)}
    if MAPS.exists():
        old = readj(MAPS / 'export_report.json')
        require(old['production_sha256'] == rh and old['spatial_source_sha256'] == input_hashes, 'Existing maps from different inputs; refusing overwrite')
        check_maps(MAPS, old)
        print('Reused verified export', flush=True)
        return
    melt = np.zeros((N, 14, 21))
    rain = np.zeros_like(melt)
    cover = np.zeros((14, 21))
    manifest = []
    source = np.zeros(2)
    scaled = np.zeros(2)
    for i, gid in enumerate(sorted(jobs), 1):
        job = jobs[gid]
        ov = overlap[gid]
        route = route_candidate(job, overlap)
        require(route in cells and all(((r, c) in cells for _, r, c in ov)), f'Outside basin {gid}')
        area = math.fsum((a for a, _, _ in ov)) / 1000000.0
        reported = float(job['polygon_overlap_besham_km2'])
        require(abs(area - reported) <= max(3e-05, 2e-06 * max(area, 1)), f'Overlap mismatch {gid}')
        if gid not in ok:
            continue
        d = PROD / gid
        r = p.verify_existing(d, job, prov)
        values = read_daily(d / 'daily.csv')
        tot = values.sum(axis=0)
        e = entries[gid]
        require(close(tot, [e['melt_m3'], e['rain_m3']]) and close(tot, [r['diagnostic_summary']['total_disaggregated_melt_m3'], r['diagnostic_summary']['total_disaggregated_rain_m3']]), f'Volume mismatch {gid}')
        fraction = add_glacier(melt, rain, cover, job, ov, route, values)
        source += tot
        scaled += tot * fraction
        manifest.append(dict(rgi_id=gid, release_row=route[0], release_col=route[1], release_verified=False, release_class=job['release_status'], overlap_fraction=fraction, overlap_km2=area, daily_sha256=r['daily_csv_sha256'], scaled_melt_m3=float(tot[0] * fraction), scaled_rain_m3=float(tot[1] * fraction)))
        if i % 500 == 0:
            print(f'Aggregated {i}/{len(jobs)}', flush=True)
    require(close(source, [report['total_melt_m3'], report['total_rain_m3']]) and close(scaled, [melt.sum(), rain.sum()]), 'Basin volume closure')
    cf = np.zeros_like(cover)
    for (r, c), (lat, lon, cellarea) in cells.items():
        y, x = regional_index(r, c)
        require(np.isclose(lat, LAT[y]) and np.isclose(lon, LON[x]), 'Basin grid mismatch')
        require(cover[y, x] <= cellarea * (1 + 1e-05), 'Cover exceeds cell')
        cf[y, x] = min(1.0, cover[y, x] / cellarea)
    require(np.isclose(cover.sum() / 1000000.0, area_total - missing_area, rtol=1e-06, atol=0.05), 'Cover overlap closure')
    MAPS.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.era5-export-', dir=MAPS.parent))
    try:
        write_regional(stage / NAMES['melt'], VARIABLES['melt'], melt, meteo)
        write_regional(stage / NAMES['rain'], VARIABLES['rain'], rain, meteo)
        write_netcdf(stage / NAMES['cover'], cf, LAT, LON, name='on_area', units='1', yearly=True)
        with (stage / 'manifest.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(manifest[0]))
            w.writeheader()
            w.writerows(manifest)
        result = dict(status='ERA5_REGIONAL_GLACIER_FORCING_PROVISIONAL', production_sha256=rh, spatial_source_sha256=input_hashes, glaciers=len(ok), omitted_glaciers=report['blocked_entries'], missing_overlap_km2=missing_area, coverage_percent=100 * (1 - missing_area / area_total), represented_polygon_overlap_km2=float(cover.sum() / 1000000.0), source_m3=dict(zip(('melt', 'rain'), map(float, source))), exported_m3=dict(zip(('melt', 'rain'), map(float, scaled))), output_sha256={name: p.sha(stage / name) for name in NAMES.values()}, manifest_sha256=p.sha(stage / 'manifest.csv'), routing='Largest in-mask polygon overlap, not verified terminus; capped area scaling; static cover for included glaciers only')
        check_maps(stage, result)
        writej(stage / 'export_report.json', result)
        require(p.sha(report_path) == rh, 'Production report changed during export')
        require(not MAPS.exists(), 'Output appeared during export')
        stage.rename(MAPS)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    print('EXPORT COMPLETE', MAPS, flush=True)

def self_test():
    m = np.zeros((2, 14, 21))
    r = np.zeros_like(m)
    c = np.zeros((14, 21))
    job = {'rgi_id': 'test', 'polygon_overlap_besham_km2': '1', 'rgi_area_km2': '2'}
    f = add_glacier(m, r, c, job, [(1000000.0, 106, 505)], (106, 505), np.array([[10.0, 2.0], [20.0, 4.0]]))
    assert f == 0.5 and m.sum() == 15 and (r.sum() == 3) and (c.sum() == 1000000.0)
    try:
        regional_index(105, 504)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Outside cut accepted')
    print('EXPORT SELF TEST PASS')

from collections import defaultdict
from datetime import date
OVERLAP=BASIN/'overlap.csv'
CELLS=BASIN/'cells.csv'
ROUTES=BASIN/'routing.csv'
WORKLIST=BASIN/'worklist.csv'
class BridgeError(RuntimeError): pass
class FinalError(RuntimeError): pass

def csvread(path):
    return read_table(path)

def read_csv(path):
    if not path.is_file():
        raise BridgeError(f'Missing required source file: {path}')
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise BridgeError(f'Missing CSV header: {path}')
        return [dict(r) for r in reader], set(reader.fieldnames)

def as_float(value, label, nonnegative=True):
    try:
        x = float(value)
    except (TypeError, ValueError):
        raise BridgeError(f'Invalid {label}: {value!r}') from None
    if not math.isfinite(x) or (nonnegative and x < 0):
        raise BridgeError(f'Nonfinite/negative {label}: {value!r}')
    return x

def load_grid(root):
    """Construct a compact, EXACT 0.5-degree CWatM-aligned bbox, using recorded cell areas."""
    basin, bfields = read_csv(root / CELLS)
    route, rfields = read_csv(root / ROUTES)
    req_b = {'row', 'col', 'latitude', 'longitude', 'cell_area_km2'}
    req_r = {'grid_row', 'grid_col', 'routes_to_provisional_besham_outlet'}
    if not req_b <= bfields or not req_r <= rfields:
        raise BridgeError('Unexpected basin-cell or routing CSV schema')
    if len(basin) != 67 or len(route) != 67:
        raise BridgeError('Expected 67 provisional Besham CWatM cells')
    cells = {}
    for r in basin:
        key = (int(r['row']), int(r['col']))
        if key in cells:
            raise BridgeError(f'Duplicate basin grid cell: {key}')
        area = as_float(r['cell_area_km2'], 'cell area') * 1e6
        if not area:
            raise BridgeError('Zero basin cell area')
        cells[key] = (as_float(r['latitude'], 'latitude', False),
                      as_float(r['longitude'], 'longitude', False), area)
    dest = set()
    for r in route:
        key = int(r['grid_row']), int(r['grid_col'])
        if key in dest or r['routes_to_provisional_besham_outlet'].strip().lower() not in ('true', '1'):
            raise BridgeError('Duplicate or unconnected CWatM routing cell')
        dest.add(key)
    if dest != set(cells):
        raise BridgeError('Basin cell set differs from verified routing CSV')
    row_min, row_max = min(r for r, _ in cells), max(r for r, _ in cells)
    col_min, col_max = min(c for _, c in cells), max(c for _, c in cells)
    nlat, nlon = row_max - row_min + 1, col_max - col_min + 1
    # Global latitude/longitude cell centers, inferred from the observed 0.5° grid.
    anchor_key = sorted(cells)[0]
    lat_anchor, lon_anchor, _ = cells[anchor_key]
    lats = [lat_anchor - (r - anchor_key[0]) * .5 for r in range(row_min, row_max+1)]
    lons = [lon_anchor + (c - anchor_key[1]) * .5 for c in range(col_min, col_max+1)]
    for (r, c), (lat, lon, _) in cells.items():
        if abs(lat - lats[r-row_min]) > 1e-5 or abs(lon-lons[c-col_min]) > 1e-5:
            raise BridgeError('Grid is not an aligned descending-latitude 0.5-degree grid')
    return cells, (row_min, col_min), lats, lons

def write_netcdf(path, data, lat, lon, *, name, units, yearly=False):
    """Write a complete 360x720 CWatM-aligned NetCDF4 with compressed time slices.

    Global cover grid is required by CWatM land-cover reads (cut=False).
    The daily melt and rain maps use the separate regional writer.
    """
    import numpy as np
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise BridgeError('netCDF4 package unavailable in the active OGGM environment') from exc
    latg = 89.75 - np.arange(360, dtype='f8') * .5
    long = -179.75 + np.arange(720, dtype='f8') * .5
    # The compact basin array must match a contiguous cell-aligned portion of
    # the native global CWatM 30-minute grid, never silently reprojected.
    lat_idx = [round((89.75 - x)*2) for x in lat]
    lon_idx = [round((x + 179.75)*2) for x in lon]
    if (min(lat_idx)<0 or max(lat_idx)>=360 or min(lon_idx)<0 or max(lon_idx)>=720
        or any(abs(latg[i]-x)>1e-6 for i,x in zip(lat_idx,lat))
        or any(abs(long[i]-x)>1e-6 for i,x in zip(lon_idx,lon))
        or lat_idx != list(range(lat_idx[0],lat_idx[0]+len(lat)))
        or lon_idx != list(range(lon_idx[0],lon_idx[0]+len(lon)))):
        raise BridgeError('Source grid does not align with the native global 0.5-degree grid')
    r0, r1 = lat_idx[0], lat_idx[-1] + 1
    c0, c1 = lon_idx[0], lon_idx[-1] + 1
    ntime = 15 if yearly else 5478
    with Dataset(str(path), 'w', format='NETCDF4_CLASSIC') as nc:
        nc.createDimension('time', ntime)
        nc.createDimension('lat', 360)
        nc.createDimension('lon', 720)
        t=nc.createVariable('time','f8',('time',))
        t.units='days since 2005-01-01 00:00:00'
        t.calendar='standard'
        t.standard_name='time'
        if yearly:
            t[:]=np.array([(date(y,1,1)-date(2005,1,1)).days for y in range(2005,2020)],dtype='f8')
        else:
            t[:]=np.arange(5478,dtype='f8')
        yl=nc.createVariable('lat','f8',('lat',))
        xl=nc.createVariable('lon','f8',('lon',))
        yl.units='degrees_north'; yl.standard_name='latitude'; yl[:]=latg
        xl.units='degrees_east'; xl.standard_name='longitude'; xl[:]=long
        v=nc.createVariable(name,'f4',('time','lat','lon'),zlib=True,complevel=3,
                            shuffle=True,chunksizes=(1,360,720))
        v.units=units
        v.long_name={'melt_on_glacier':'Daily glacier melt volume at assumed release cells',
                     'liq_prcp_on_glacier':'Daily glacier rainfall volume at assumed release cells',
                     'on_area':'Fraction of each CWatM grid cell covered by included glaciers'}[name]
        nc.title='PROVISIONAL OGGM-CWatM coupling input - NOT basin-complete glacier runoff'
        nc.source='Custom daily OGGM-anchored surface model'
        nc.approved_for_flood_analysis='NO'
        # One compressed global frame at a time. No global (time,lat,lon)
        # allocation, so even a 15-year daily dataset uses little RAM.
        frame=np.zeros((360,720),dtype='f4')
        if yearly:
            frame[r0:r1,c0:c1]=np.asarray(data,dtype='f4')
            for i in range(ntime):
                v[i,:,:]=frame
        else:
            for i in range(ntime):
                frame[r0:r1,c0:c1]=np.asarray(data[i,:,:],dtype='f4')
                v[i,:,:]=frame

def load_overlap(root):
    overlap = csvread(root / OVERLAP)
    by = defaultdict(list)
    for row in overlap:
        try:
            area = float(row['glacier_polygon_area_in_cell_m2'])
            rr = int(row['grid_row'])
            cc = int(row['grid_col'])
        except (KeyError, ValueError) as e:
            raise FinalError('Invalid polygon-grid overlap CSV') from e
        if area <= 0:
            raise FinalError('Nonpositive polygon overlap')
        by[row['rgi_id']].append((area, rr, cc))
    return by

def route_candidate(job,overlap):
 gid=job['rgi_id']
 if gid not in overlap:raise FinalError(f'Glacier has no polygon overlap: {gid}')
 area,row,col=max(overlap[gid],key=lambda x:x[0])
 if job['release_status']=='PROVISIONAL_ONE_CELL_UNVERIFIED':
  if len(overlap[gid])!=1:raise FinalError(f'Original single-cell classification disagrees with overlap: {gid}')
  if (row,col)!=(int(job['release_grid_row']),int(job['release_grid_col'])):raise FinalError(f'Original single-cell row/column mismatch: {gid}')
 return row,col


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--accept-omissions', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    with (PROD / '.production.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        export(args.accept_omissions)

if __name__=='__main__':main()
