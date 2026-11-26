#!/usr/bin/env python3
"""Resumable acquisition of missing Upper Indus ERA5 hours from ECMWF ARCO.
Run from project root. Raw hourly data only; no model changes.
"""
import argparse, calendar, hashlib, json, os, time
from pathlib import Path
import numpy as np

from .common import ROOT
CDS = ROOT/'data/climate/era5/raw/reanalysis-era5-single-levels'
OUT = ROOT/'data/climate/era5/raw/ecmwf_arco'
URL = ('https://arco.datastores.ecmwf.int/cadl-arco-geo-002/arco/'
       'reanalysis_era5_single_levels/sfc/geoChunked.zarr')
VARS = dict(t2m='2m_temperature', d2m='2m_dewpoint_temperature',
            sp='surface_pressure', u10='10m_u_component_of_wind',
            v10='10m_v_component_of_wind', tp='total_precipitation',
            ssrd='surface_solar_radiation_downwards',
            strd='surface_thermal_radiation_downwards')
UNITS = dict(t2m='K',d2m='K',sp='Pa',u10='m s-1',v10='m s-1',tp='m',ssrd='J m-2',strd='J m-2')
AREA = [38.,71.5,30.,83.]
START = np.datetime64('2005-01-01T00','h')
END = np.datetime64('2020-01-01T00','h')
TIMES = np.arange(START, END+np.timedelta64(1,'h'), np.timedelta64(1,'h'))
LAT = np.arange(38.,29.99,-.25)
LON = np.arange(71.5,83.01,.25)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def json_write(path, obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    os.replace(tmp,path)

def cds_valid(y,m,boundary=False):
    stem=f'era5_single_levels_{y:04d}_{m:02d}'+('_boundary00' if boundary else '')
    file=CDS/(stem+'.grib'); meta=CDS/(stem+'.request.json')
    if not file.exists(): return None
    try:
        info=json.loads(meta.read_text()); r=info['request']
        days=['01'] if boundary else [f'{d:02d}' for d in range(1,calendar.monthrange(y,m)[1]+1)]
        hours=['00:00'] if boundary else [f'{h:02d}:00' for h in range(24)]
        assert info['dataset']=='reanalysis-era5-single-levels'
        assert r['year']==[f'{y:04d}'] and r['month']==[f'{m:02d}']
        assert sorted(r['day'])==days and sorted(r['time'])==hours
        assert set(r['variable'])==set(VARS.values())
        assert r['area']==AREA and r['grid']==[.25,.25]
        assert info['size_bytes']==file.stat().st_size and info['sha256']==sha(file)
        with file.open('rb') as f:
            assert f.read(4)==b'GRIB'; f.seek(-4,2); assert f.read()==b'7777'
        # Validate actual message coverage, coordinates, and step types, not just header.
        from eccodes import codes_grib_new_from_file, codes_get, codes_release
        aliases={'2t':'t2m','2d':'d2m','10u':'u10','10v':'v10'}
        from datetime import datetime
        expected={int(str(np.datetime64(f'{y:04d}-{m:02d}-{d}T{h}').astype('datetime64[h]')).replace('-','').replace('T','')) for d in days for h in hours}
        seen={v:set() for v in VARS}
        with file.open('rb') as f:
            while True:
                g=codes_grib_new_from_file(f)
                if g is None: break
                try:
                    v=codes_get(g,'shortName'); v=aliases.get(v,v)
                    assert v in seen
                    date=int(codes_get(g,'validityDate')); hm=int(codes_get(g,'validityTime'))
                    assert hm%100==0
                    stamp=date*100+hm//100
                    assert stamp in expected and stamp not in seen[v]
                    seen[v].add(stamp)
                    assert (codes_get(g,'Ni'),codes_get(g,'Nj'))==(47,33)
                    assert abs(codes_get(g,'latitudeOfFirstGridPointInDegrees')-38)<1e-6
                    assert abs(codes_get(g,'longitudeOfFirstGridPointInDegrees')-71.5)<1e-6
                    assert abs(codes_get(g,'latitudeOfLastGridPointInDegrees')-30)<1e-6
                    assert abs(codes_get(g,'longitudeOfLastGridPointInDegrees')-83)<1e-6
                    assert codes_get(g,'stepType')==('accum' if v in ('tp','ssrd','strd') else 'instant')
                finally: codes_release(g)
        assert all(s==expected for s in seen.values())
        return {'file':str(file.relative_to(ROOT)), 'sha256':info['sha256'], 'year':y,'month':m,'boundary':boundary}
    except Exception as exc:
        print(f'NOT REUSED (preserved): {file.name}: {type(exc).__name__}: {exc}',flush=True)
        return None

def plan():
    missing = np.ones(len(TIMES), dtype=bool)
    reused = []
    for y in range(2005, 2021):
        for m in range(1, 13) if y < 2020 else [1]:
            a = np.datetime64(f'{y:04d}-{m:02d}', 'M')
            b = a + np.timedelta64(1, 'M')
            mask = (TIMES >= a) & (TIMES < b)
            item = cds_valid(y, m) if y < 2020 else None
            if item:
                missing[mask] = False
                reused.append(item)
                print('REUSE', item['file'], flush=True)
            else:
                item = cds_valid(y, m, True)
                if item:
                    missing[TIMES == a.astype('datetime64[h]')] = False
                    reused.append(item)
                    print('REUSE', item['file'], flush=True)
    return (missing, reused)

def runs(idx,chunk):
    # Each returned slice intersects exactly one native chunk.
    for group in np.unique(idx//chunk):
        use=idx[idx//chunk==group]
        yield slice(int(use[0]),int(use[-1])+1),use

def main():
    p=argparse.ArgumentParser(); p.add_argument('--plan',action='store_true'); args=p.parse_args()
    if Path.cwd() != ROOT: raise SystemExit('Run from project root.')
    # Ensure ecCodes works before deciding which existing files to skip.
    from eccodes import codes_grib_new_from_file
    if (ROOT/'data/climate/era5/raw').is_symlink():
        raise SystemExit('Raw ERA5 inputs are linked from another project; do not write to them. Run prepare-climate with the existing COMPLETE.json.')
    OUT.mkdir(parents=True,exist_ok=True)
    lock=OUT/'.acquisition.lock'
    import fcntl
    with lock.open('a') as lf:
        try: fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise SystemExit('Another ARCO acquisition is running.')
        run(args)

def run(args):
    missing,reused=plan()
    hours=TIMES[missing]
    spec={'schema':1,'source':URL,'area':AREA,'grid':.25,'start':str(START),'end_inclusive':str(END),'variables':list(VARS),'reused_cds':reused,'missing_hours_sha256':hashlib.sha256(hours.astype('int64').tobytes()).hexdigest(),'missing_hours':len(hours)}
    print(f'Total hours: {len(TIMES)}; reused: {len(TIMES)-len(hours)}; ARCO: {len(hours)}',flush=True)
    print(f'Uncompressed remaining float32 payload: {len(hours)*len(LAT)*len(LON)*len(VARS)*4/2**30:.2f} GiB',flush=True)
    if args.plan: return
    # Distinct plan directories preserve old partial acquisitions if CDS inventory changes.
    ident=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:16]
    dest=OUT/ident; dest.mkdir(exist_ok=True)
    json_write(dest/'acquisition_plan.json',spec)
    if not len(hours):
        json_write(dest/'COMPLETE.json',{'plan':ident,'tiles':[],'reused_cds':reused,'hours':0,'spatial_points':len(LAT)*len(LON),'variables':list(VARS)})
        print('All required hours already available in CDS.'); return
    import xarray as xr
    key=next(line.split(':',1)[1].strip() for line in Path(os.environ.get('CDSAPI_RC',str(Path.home()/'.cdsapirc'))).read_text().splitlines() if line.startswith('key:'))
    ds=xr.open_zarr(URL,consolidated=True,chunks=None,storage_options={'headers':{'Authorization':f'Bearer {key}'}})
    assert set(VARS)<=set(ds.data_vars), 'Missing ARCO variables'
    t=np.asarray(ds.time.values).astype('datetime64[h]'); yi=np.array([np.where(ds.latitude.values==v)[0].item() for v in LAT]); xi=np.array([np.where(ds.longitude.values==v)[0].item() for v in LON])
    ti=np.searchsorted(t,hours)
    assert np.all(ti<len(t)) and np.array_equal(t[ti],hours), 'ARCO missing time coverage'
    assert np.all(np.diff(ti)>0)
    # Native spatial chunks differ by variable in principle; use each independently.
    tasks=[]
    for v in VARS:
        a=ds[v]; chunks=dict(zip(a.dims,a.encoding['chunks']))
        assert set(a.dims)=={'time','latitude','longitude'},f'Unexpected dimensions: {v}: {a.dims}'
        print(v,'native chunks',chunks,flush=True)
        for ys,yidx in runs(np.sort(yi),chunks['latitude']):
            for xs,xidx in runs(np.sort(xi),chunks['longitude']):
                tasks.append((v,ys,xs,chunks['time']))
    completed=[]
    for n,(v,ys,xs,ct) in enumerate(tasks,1):
        name=f'{v}_y{ys.start}-{ys.stop}_x{xs.start}-{xs.stop}.nc'
        path=dest/name; meta=dest/(name+'.json')
        expected={'plan':ident,'variable':v,'y':[ys.start,ys.stop],'x':[xs.start,xs.stop],'hours':len(hours)}
        try:
            old=json.loads(meta.read_text())
            if all(old.get(k)==val for k,val in expected.items()) and path.stat().st_size==old['bytes'] and sha(path)==old['sha256']:
                print(f'[{n}/{len(tasks)}] SKIP {name}',flush=True); completed.append(old); continue
        except (OSError,ValueError,KeyError): pass
        print(f'[{n}/{len(tasks)}] Fetch {name}',flush=True)
        if path.exists(): raise SystemExit(f'Unverified existing tile {path}; inspect before moving aside and retrying.')
        a=ds[v]; ny=ys.stop-ys.start; nx=xs.stop-xs.start
        data=np.empty((len(hours),ny,nx),dtype=np.float32); offset=0
        for ts,indices in runs(ti,ct):
            for attempt in range(5):
                try:
                    block=a.isel(time=ts,latitude=ys,longitude=xs).transpose('time','latitude','longitude').values
                    selected=block[indices-ts.start,:,:]
                    if not np.isfinite(selected).all(): raise ValueError('Non-finite ARCO data')
                    data[offset:offset+len(indices)]=selected; offset+=len(indices)
                    break
                except Exception as e:
                    if attempt==4: raise
                    print(f'  Retry {attempt+1}/4: {type(e).__name__}',flush=True); time.sleep(5*2**attempt)
        assert offset==len(hours)
        tile=xr.Dataset({v:(('time','latitude','longitude'),data,{'units':a.attrs.get('units',UNITS[v])})},coords={'time':hours.astype('datetime64[ns]'),'latitude':ds.latitude.values[ys],'longitude':ds.longitude.values[xs]},attrs={'source':URL,'acquisition_plan':ident,'time_convention':'UTC validity time; tp, ssrd, strd are preceding-hour accumulations'})
        tmp=dest/(name+'.partial')
        tile.to_netcdf(tmp,engine='netcdf4',encoding={v:{'zlib':True,'complevel':1,'dtype':'float32'}})
        # Reopen local output before promotion; all entries checked during retrieval.
        with xr.open_dataset(tmp,engine='netcdf4') as check:
            assert check[v].shape==data.shape and np.array_equal(check.time.values.astype('datetime64[h]'),hours)
        info={**expected,'file':name,'bytes':tmp.stat().st_size,'sha256':sha(tmp)}
        # Promote metadata first; an interrupted promotion leaves no falsely accepted tile.
        json_write(meta,info); os.replace(tmp,path); completed.append(info)
        print(f'  Saved {path.stat().st_size/2**20:.2f} MiB',flush=True)
    json_write(dest/'COMPLETE.json',{'plan':ident,'reused_cds':reused,'tiles':completed,'hours':len(hours),'spatial_points':len(LAT)*len(LON),'variables':list(VARS)})
    print('ACQUISITION COMPLETE:',dest/'COMPLETE.json',flush=True)
    ds.close()

if __name__=='__main__': main()
