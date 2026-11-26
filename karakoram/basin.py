"""Delineate the fixed Besham outlet and intersect RGI62 glacier polygons.

The outlet is a documented research choice, not an automatic gauge fit. The
polygon intersection is an area proxy; it does not establish glacier drainage.
"""
from __future__ import annotations
import argparse, csv, math, re, os, shutil, tempfile
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from .common import ROOT, CONFIG, SPATIAL, BASIN, OUTLET, require, sha, read_csv, write_csv, atomic_json, read_json

CODE_TO_SOUTH_EAST = {1:(1,-1), 2:(1,0), 3:(1,1), 4:(0,-1), 6:(0,1), 7:(-1,-1), 8:(-1,0), 9:(-1,1)}
DIRS = CODE_TO_SOUTH_EAST
RGI_ID = re.compile(r'^RGI60-(\d{2})\.\d{5}$')
class InventoryError(RuntimeError): pass
class RoutingScreenError(RuntimeError): pass

@dataclass
class Network:
    downstream: np.ndarray  # flat downstream index, -1 for pits/boundaries
    area_m2: np.ndarray      # flat cell area
    upstream_area_m2: np.ndarray
    valid: np.ndarray
    shape: tuple[int, int]
    total_nodes: int
    sinks: int
    terminated_at_border_or_missing: int
    ordered_nodes: int

def build_network(ldd: np.ndarray, cellarea: np.ndarray, lat: np.ndarray,
                  lon: np.ndarray) -> Network:
    """Trace LDD's one-way graph and accumulate upstream area, including outlet cells.

    Rejects any network cycles; does not silently repair topology or write maps.
    """
    if ldd.shape != cellarea.shape or ldd.shape != (len(lat), len(lon)):
        raise ValueError("Grid shapes do not agree")
    if np.any(np.isfinite(ldd) & (ldd != np.rint(ldd))):
        raise ValueError("LDD has non-integer direction codes")

    codes = np.where(np.isfinite(ldd), ldd, 0).astype(np.int16)
    invalid_codes = set(np.unique(codes)) - set(range(10))
    if invalid_codes:
        raise ValueError(f"Invalid LDD direction codes: {sorted(invalid_codes)}")
    # Codes 1..9 are valid where local cell area exists; code 0 is no data.
    valid_2d = (codes >= 1) & (codes <= 9) & np.isfinite(cellarea) & (cellarea > 0)
    valid = valid_2d.ravel()
    rows, cols = codes.shape
    n = rows * cols
    dr_south = 1 if lat[0] > lat[-1] else -1
    dc_east = 1 if lon[-1] > lon[0] else -1
    downstream = np.full(n, -1, dtype=np.int64)
    border_ends = 0

    for code, (south, east) in CODE_TO_SOUTH_EAST.items():
        rr, cc = np.where((codes == code) & valid_2d)
        r2 = rr + dr_south * south
        c2 = cc + dc_east * east
        bounds = (r2 >= 0) & (r2 < rows) & (c2 >= 0) & (c2 < cols)
        target_valid = np.zeros(rr.size, dtype=bool)
        target_valid[bounds] = valid_2d[r2[bounds], c2[bounds]]
        good = bounds & target_valid
        downstream[rr[good] * cols + cc[good]] = r2[good] * cols + c2[good]
        border_ends += int((~good).sum())

    indegree = np.bincount(downstream[downstream >= 0], minlength=n).astype(np.int64)
    queue = deque(int(i) for i in np.flatnonzero(valid & (indegree == 0)))
    area_flat = np.where(valid, cellarea.ravel(), 0.0).astype(np.float64)
    accumulated = area_flat.copy()
    processed = 0
    while queue:
        cell = queue.popleft()
        processed += 1
        next_cell = int(downstream[cell])
        if next_cell >= 0:
            accumulated[next_cell] += accumulated[cell]
            indegree[next_cell] -= 1
            if indegree[next_cell] == 0:
                queue.append(next_cell)
    total_nodes = int(valid.sum())
    if processed != total_nodes:
        raise ValueError(
            f"LDD graph contains cycles or unresolved nodes: {total_nodes - processed} / {total_nodes}. "
            "Do not delineate without checking the source routing map."
        )
    return Network(downstream, area_flat, accumulated, valid, (rows, cols),
                   total_nodes, int(((codes == 5) & valid_2d).sum()), border_ends, processed)

def upstream_mask(downstream: np.ndarray, valid: np.ndarray, outlet: int) -> np.ndarray:
    """Traverse upstream edges of an acyclic LDD network using O(N) storage/time.

    The contributing catchment *includes* the outlet grid cell. The source
    network's topological validation is performed by build_network().
    """
    n = len(downstream)
    if valid.shape != (n,) or not 0 <= outlet < n or not bool(valid[outlet]):
        raise ValueError('Selected outlet is not a valid LDD cell')
    parents = np.flatnonzero(downstream >= 0)
    if np.any(downstream[parents] >= n):
        raise ValueError('Out-of-range downstream index')
    first_child = np.full(n, -1, dtype=np.int64)
    next_sibling = np.full(n, -1, dtype=np.int64)
    for child in parents:
        receiver = int(downstream[child])
        next_sibling[child] = first_child[receiver]
        first_child[receiver] = int(child)
    mask = np.zeros(n, dtype=bool)
    stack = [int(outlet)]
    while stack:
        cell = stack.pop()
        if mask[cell]:
            raise ValueError('Unexpected repeated node/cycle in upstream traversal')
        mask[cell] = True
        child = int(first_child[cell])
        while child != -1:
            stack.append(child)
            child = int(next_sibling[child])
    return mask

def make_cell_polygons(mask: np.ndarray, lat: np.ndarray, lon: np.ndarray):
    """Build actual CWatM geographic cell polygons; no centroid-only assignment."""
    from shapely.geometry import box
    result = []
    for r, c in np.argwhere(mask):
        r, c = int(r), int(c)
        geom = box(float(lon[c]) - 0.25, float(lat[r]) - 0.25,
                   float(lon[c]) + 0.25, float(lat[r]) + 0.25)
        result.append((r, c, geom))
    return result

def examine_glaciers(gdf, region: str, cells: list,
                     *, transformer=None) -> tuple[list[dict], list[dict], dict]:
    """Intersect glacier polygons with active cells, retaining uncertain edges.

    Areas use a projected equal-area CRS (EPSG:6933) for *spatial fractions*
    only. Polygon overlap is NOT proof of drainage outlet or water routing.
    """
    from shapely.ops import transform
    from shapely.ops import unary_union
    from pyproj import Transformer
    if transformer is None:
        transformer = Transformer.from_crs('EPSG:4326', 'EPSG:6933', always_xy=True).transform
    cell_entries = [(r, c, poly.bounds, transform(transformer, poly)) for r, c, poly in cells]
    basin = unary_union([row[3] for row in cell_entries])
    selected = []
    allocation = []
    checked = 0
    for _, rec in gdf.iterrows():
        geom = rec.geometry
        if geom is None or geom.is_empty:
            raise InventoryError('RGI contains an empty glacier geometry')
        if not geom.is_valid:
            raise InventoryError(f'Invalid glacier geometry; refusing silent repair: {rec.get("RGIId")}')
        raw_id = str(rec.get('RGIId', ''))
        match = RGI_ID.fullmatch(raw_id)
        if not match or match.group(1) != region:
            raise InventoryError(f'RGI62 source has incompatible ID/region: {raw_id}')
        try:
            rgi_area = float(rec.get('Area'))
        except (ValueError, TypeError):
            raise InventoryError(f'Invalid RGI area for {raw_id}') from None
        if not math.isfinite(rgi_area) or rgi_area <= 0:
            raise InventoryError(f'RGI glacier {raw_id} has nonpositive/nonfinite area')
        checked += 1
        geom_eq = transform(transformer, geom)
        a_full = float(geom_eq.area)
        if not math.isfinite(a_full) or a_full <= 0:
            raise InventoryError(f'Invalid projected area for {raw_id}')
        a_basin = float(geom_eq.intersection(basin).area)
        if a_basin <= 1e-3:
            continue  # Touching the basin boundary at a line/point is not area overlap.
        frac = min(1.0, max(0.0, a_basin / a_full))
        # Guard against tiny floating-point perturbations at boundaries.
        if frac > 1 + 1e-6:
            raise InventoryError(f'Overlap fraction exceeds 1 for {raw_id}')
        bbox = geom.bounds
        cell_count = 0
        sum_cells = 0.0
        for r, c, cbbox, cpoly in cell_entries:
            if bbox[2] <= cbbox[0] or bbox[0] >= cbbox[2] or bbox[3] <= cbbox[1] or bbox[1] >= cbbox[3]:
                continue
            overlap = float(geom_eq.intersection(cpoly).area)
            if overlap <= 1e-3:
                continue
            allocation.append({'rgi_id': raw_id, 'grid_row': r, 'grid_col': c,
                               'glacier_polygon_area_in_cell_m2': f'{overlap:.3f}',
                               'fraction_of_glacier_polygon': f'{overlap / a_full:.10f}'})
            sum_cells += overlap
            cell_count += 1
        if cell_count == 0 or not math.isclose(sum_cells, a_basin, abs_tol=1.0, rel_tol=1e-7):
            raise InventoryError(f'Cell overlap does not close for {raw_id}: {sum_cells}, {a_basin}')
        selected.append({'rgi_id': raw_id, 'rgi_region': region,
                         'rgi_area_km2': f'{rgi_area:.6f}',
                         'polygon_area_km2_equal_area': f'{a_full / 1e6:.6f}',
                         'polygon_overlap_besham_km2_equal_area': f'{a_basin / 1e6:.6f}',
                         'fraction_polygon_inside_besham': f'{frac:.10f}',
                         'rgi_area_inside_besham_km2_estimate': f'{rgi_area * frac:.6f}',
                         'membership': 'within_mask' if frac >= (1 - 1e-7) else 'boundary_crossing',
                         'overlap_cell_count': cell_count,
                         'rgi_centre_lat': str(rec.get('CenLat', '')),
                         'rgi_centre_lon': str(rec.get('CenLon', ''))})
    return selected, allocation, {'screened_bbox_glaciers': checked,
                                  'included_glaciers': len(selected),
                                  'boundary_crossing_glaciers': sum(x['membership'] == 'boundary_crossing' for x in selected)}

def read_rgi(path: Path, basin_bbox: tuple[float, float, float, float]):
    """Load RGI62 geometries only in basin bounding rectangle, with CRS check."""
    import geopandas as gpd
    from pyproj import CRS
    # GeoPandas supports pyogrio and Fiona readers. Use pyogrio explicitly:
    # the project's OGGM environment can write shapefiles using pyogrio,
    # but need not have the optional Fiona package installed.
    gdf = gpd.read_file(path, bbox=basin_bbox, engine='pyogrio')
    if gdf.crs is None:
        raise InventoryError('RGI read yielded geometry with no coordinate reference system')
    if not CRS.from_user_input(gdf.crs).equals(CRS.from_epsg(4326)):
        raise InventoryError('RGI feature CRS does not match EPSG:4326')
    if not {'RGIId', 'Area'}.issubset(gdf.columns):
        raise InventoryError('RGI shapefile missing RGIId or Area fields')
    if gdf['RGIId'].duplicated().any():
        raise InventoryError('RGI shapefile contains duplicate glacier IDs')
    return gdf

def finite_number(value: str, what: str, *, positive: bool=False) -> float:
    try:
        number=float(value)
    except (TypeError,ValueError) as exc:
        raise RoutingScreenError(f'Invalid {what}: {value!r}') from exc
    if not math.isfinite(number) or (positive and number <= 0):
        raise RoutingScreenError(f'Invalid {what}: {value!r}')
    return number

def cell_index(value: str, bound: int, what: str) -> int:
    try:
        n=int(value)
    except (TypeError,ValueError) as exc:
        raise RoutingScreenError(f'Invalid {what}: {value!r}') from exc
    if str(n) != str(value).strip() or not 0 <= n < bound:
        raise RoutingScreenError(f'Invalid {what}: {value!r}')
    return n

def trace_besham_mask(mask: np.ndarray, ldd: np.ndarray) -> dict[tuple[int,int], dict]:
    """Prove each of the 67 mask cells reaches the outlet through in-mask LDD.

    Uses the original PCRaster keypad directions (latitude decreases by row),
    not glacier polygon overlap or an inferred glacier terminus. This is a
    separate 67-cell trace, not a full native CWatM hydrological run.
    """
    if mask.shape != ldd.shape or not mask[OUTLET] or int(mask.sum()) != 67:
        raise RoutingScreenError('Expected 67-cell Besham mask and aligned LDD')
    results={}
    nrows,ncols=mask.shape
    for rr,cc in np.argwhere(mask):
        start=(int(rr),int(cc))
        seen=set()
        cursor=start
        moves=0
        while cursor != OUTLET:
            if cursor in seen:
                raise RoutingScreenError(f'CWatM LDD loop from mask cell {start}')
            seen.add(cursor)
            r,c=cursor
            code=ldd[r,c]
            if not math.isfinite(float(code)) or code != int(code) or int(code) not in DIRS:
                raise RoutingScreenError(f'Invalid LDD direction {code} at {cursor}')
            dr,dc=DIRS[int(code)]
            cursor=(r+dr,c+dc)
            if not (0<=cursor[0]<nrows and 0<=cursor[1]<ncols and mask[cursor]):
                raise RoutingScreenError(f'Mask cell {start} routes outside Besham mask via {cursor}')
            moves+=1
            if moves>=int(mask.sum()):
                raise RoutingScreenError(f'LDD not convergent from {start}')
        code=ldd[start]
        if not math.isfinite(float(code)) or code!=int(code) or not 1<=int(code)<=9:
            raise RoutingScreenError(f'Invalid outlet/cell LDD code for {start}')
        dr,dc=DIRS.get(int(code),(0,0))
        downstream=(start[0]+dr,start[1]+dc) if start!=OUTLET else (None,None)
        results[start]={'grid_row':start[0],'grid_col':start[1],
             'ldd_code':int(code),'downstream_row':downstream[0],
             'downstream_col':downstream[1], 'steps_to_besham_outlet':moves,
             'routes_to_provisional_besham_outlet':True}
    return results

def verify_candidate_tables(candidates: list[dict], overlap: list[dict],
                            mask: np.ndarray) -> tuple[list[dict], dict]:
    """Check every glacier/grid area allocation against polygon intersection output.

    Deliberately never converts polygon fractional overlap to melt fraction.
    Whole-glacier runoff can only be assigned after a release point is known.
    """
    cand={}
    for row in candidates:
        gid=row['rgi_id']
        if not gid.startswith('RGI60-') or gid in cand:
            raise RoutingScreenError(f'Invalid or repeated glacier ID: {gid}')
        if row['membership'] not in ('within_mask','boundary_crossing'):
            raise RoutingScreenError(f'Unexpected membership for {gid}')
        frac=finite_number(row['fraction_polygon_inside_besham'], 'glacier fraction')
        area=finite_number(row['polygon_overlap_besham_km2_equal_area'],'overlap area')
        full=finite_number(row['polygon_area_km2_equal_area'],'full polygon area',positive=True)
        rgi=finite_number(row['rgi_area_km2'],'RGI area',positive=True)
        estimated=finite_number(row['rgi_area_inside_besham_km2_estimate'],'estimated RGI area')
        # Polygon intersection rounds km² to 6 decimals and fractions to 10. Tiny
        # boundary slivers can round to zero without being invalid.
        if not (0 <= frac <= 1+1e-9) or area < 0 or estimated < 0:
            raise RoutingScreenError(f'Invalid polygon overlap for {gid}')
        area_rounding_tol=0.51e-6  # km², the exported 6-decimal precision
        if abs(area - full*frac)>area_rounding_tol*(2+abs(frac))+full*2e-9:
            raise RoutingScreenError(f'Inconsistent polygon overlap for {gid}')
        if abs(estimated-rgi*frac)>5e-5:
            raise RoutingScreenError(f'Inconsistent reported RGI overlap area for {gid}')
        if row['membership']=='within_mask' and frac<1-1e-7:
            raise RoutingScreenError(f'Interior glacier fraction not unity for {gid}')
        if row['membership']=='boundary_crossing' and frac>=1-1e-7:
            raise RoutingScreenError(f'Boundary glacier fraction not smaller than unity for {gid}')
        n=cell_index(row['overlap_cell_count'],100,'overlap_cell_count')
        if n<1:
            raise RoutingScreenError(f'Empty glacier overlap {gid}')
        cand[gid]=(row,frac,area,n)
    groups=defaultdict(list)
    seen_pairs=set()
    for row in overlap:
        gid=row['rgi_id']
        if gid not in cand:
            raise RoutingScreenError(f'Cell overlap has unknown glacier: {gid}')
        r=cell_index(row['grid_row'],mask.shape[0],'grid_row')
        c=cell_index(row['grid_col'],mask.shape[1],'grid_col')
        if not mask[r,c] or (gid,r,c) in seen_pairs:
            raise RoutingScreenError(f'Outside mask / repeated glacier cell: {gid} {r},{c}')
        seen_pairs.add((gid,r,c))
        area=finite_number(row['glacier_polygon_area_in_cell_m2'],'glacier cell area',positive=True)
        fraction=finite_number(row['fraction_of_glacier_polygon'],'glacier cell fraction')
        if not 0<=fraction<=1+1e-9:
            raise RoutingScreenError(f'Glacier cell fraction over 1: {gid}')
        groups[gid].append((r,c,area,fraction))
    output=[]
    counts=Counter()
    for gid,(row,frac,total_km2,n) in sorted(cand.items()):
        values=groups.get(gid,[])
        if len(values)!=n:
            raise RoutingScreenError(f'Glacier overlap cell count mismatch for {gid}')
        # Areas stored to 0.001 m²; total overlap to 0.000001 km².
        if abs(sum(v[2] for v in values)-total_km2*1e6)>max(1.0,1.0e-6*total_km2*1e6):
            raise RoutingScreenError(f'Glacier cell-area sum mismatch for {gid}')
        if abs(sum(v[3] for v in values)-frac)>2e-6:
            raise RoutingScreenError(f'Glacier cell-fraction sum mismatch for {gid}')
        largest=max(values,key=lambda x:(x[2],-x[0],-x[1]))
        if row['membership']=='boundary_crossing':
            cls='boundary_review_required'
        elif len(values)==1:
            cls='single_cell_candidate_release_unverified'
        else:
            cls='multi_cell_release_point_needed'
        counts[cls]+=1
        output.append({'rgi_id':gid, 'membership':row['membership'],
          'overlap_cell_count':len(values),
          'polygon_fraction_inside_mask':f'{frac:.10f}',
          'screen_class':cls,
          'diagnostic_largest_overlap_row':largest[0],
          'diagnostic_largest_overlap_col':largest[1],
          'diagnostic_largest_overlap_fraction':f'{largest[3]:.10f}',
          'single_cell_candidate_row':largest[0] if cls.startswith('single_cell') else '',
          'single_cell_candidate_col':largest[1] if cls.startswith('single_cell') else '',
          'glacier_release_point_verified':'NO',
          'whole_glacier_besham_membership_verified':'NO'})
    if len(output)!=len(candidates) or len(seen_pairs)!=len(overlap):
        raise RoutingScreenError('Glacier candidate/allocation completeness failed')
    return output,dict(counts)


def read_grid(path):
    """Read the single two-dimensional global routing raster and coordinates."""
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        lat = np.asarray(ds['lat'][:], dtype=float)
        lon = np.asarray(ds['lon'][:], dtype=float)
        fields = [v for v in ds.variables if ds[v].dimensions == ('lat','lon')]
        require(len(fields)==1, f'Ambiguous raster variable: {path}')
        data = np.ma.asarray(ds[fields[0]][:]).astype(float).filled(np.nan)
    require(np.array_equal(lat, 89.75-np.arange(360)*.5) and
            np.array_equal(lon, -179.75+np.arange(720)*.5), f'Unexpected global grid: {path}')
    return data, lat, lon


def make_jobs(candidates, screen):
    """Preserve the exact original job schema used by the glacier calculation."""
    by_id = {r['rgi_id']: r for r in screen}
    jobs = []
    for i, c in enumerate(sorted(candidates, key=lambda r:r['rgi_id'])):
        s = by_id[c['rgi_id']]
        single = s['screen_class']=='single_cell_candidate_release_unverified'
        jobs.append(dict(rgi_id=c['rgi_id'], batch_id=i//25, rgi_region=c['rgi_region'],
                         rgi_centre_lat=c['rgi_centre_lat'], rgi_centre_lon=c['rgi_centre_lon'],
                         rgi_area_km2=c['rgi_area_km2'], polygon_overlap_besham_km2=c['polygon_overlap_besham_km2_equal_area'],
                         membership=s['membership'], screen_class=s['screen_class'],
                         release_grid_row=s['single_cell_candidate_row'] if single else '',
                         release_grid_col=s['single_cell_candidate_col'] if single else '',
                         release_status='PROVISIONAL_ONE_CELL_UNVERIFIED' if single else 'RELEASE_AND_OR_MEMBERSHIP_REVIEW_REQUIRED',
                         glacier_release_verified='false', whole_glacier_membership_verified='false'))
    return jobs


def prepare(fetch_rgi=False):
    """Rebuild the mask, all candidate intersections and the glacier worklist."""
    from netCDF4 import Dataset
    ldd_path = SPATIAL/'routing/ldd.nc'
    area_path = SPATIAL/'routing/cellarea.nc'
    ldd, lat, lon = read_grid(ldd_path)
    area, _, _ = read_grid(area_path)
    network = build_network(ldd, area, lat, lon)
    mask = upstream_mask(network.downstream, network.valid, OUTLET[0]*720+OUTLET[1]).reshape(360,720)
    require(int(mask.sum())==CONFIG['expected_basin_cells'], 'Basin cell count differs from reference')
    routes = trace_besham_mask(mask, ldd)
    cells = [dict(row=int(r), col=int(c), latitude=f'{lat[r]:.5f}', longitude=f'{lon[c]:.5f}', cell_area_km2=f'{area[r,c]/1e6:.7f}') for r,c in np.argwhere(mask)]
    reference = read_csv(ROOT/'reference/basin_cells.csv')
    require([(int(x['row']),int(x['col'])) for x in cells]==[(int(x['row']),int(x['col'])) for x in reference], 'Catchment differs from frozen reference')
    for a,b in zip(cells,reference):
        require(abs(float(a['cell_area_km2'])-float(b['cell_area_km2']))<1e-5, 'Cell area differs from reference')
    polygons = make_cell_polygons(mask,lat,lon)
    bbox = (min(p.bounds[0] for _,_,p in polygons),min(p.bounds[1] for _,_,p in polygons),max(p.bounds[2] for _,_,p in polygons),max(p.bounds[3] for _,_,p in polygons))
    candidates, overlaps, sources = [], [], {}
    if fetch_rgi:
        import oggm
        from oggm import cfg, utils
        require(oggm.__version__==CONFIG['oggm_version'], 'Unexpected OGGM version')
        cfg.initialize(logging_level='WARNING')
        cfg.PARAMS['use_multiprocessing']=False
    for region in CONFIG['rgi_regions']:
        path = None
        for root in (ROOT/'data/raw/rgi', ROOT/'.runtime-home/OGGM/rgi'):
            matches=sorted(root.rglob(f'{region}_rgi62_*.shp'))
            if len(matches)==1:
                path=matches[0];break
            require(len(matches)<=1, f'Ambiguous RGI {region} source')
        if path is None and fetch_rgi:
            path=Path(utils.get_rgi_region_file(region,version='62'))
        require(path is not None and path.is_file(), f'RGI62 region {region} missing; use basin --fetch-rgi')
        sources[region]={ext:sha(path.with_suffix(ext)) for ext in ('.shp','.shx','.dbf','.prj')}
        gdf=read_rgi(path,bbox)
        c,o,info=examine_glaciers(gdf,region,polygons)
        candidates.extend(c);overlaps.extend(o)
        print('RGI',region,info,flush=True)
    candidates.sort(key=lambda r:r['rgi_id'])
    overlaps.sort(key=lambda r:(r['rgi_id'],int(r['grid_row']),int(r['grid_col'])))
    screen,counts=verify_candidate_tables(candidates,overlaps,mask)
    jobs=make_jobs(candidates,screen)
    require(len(jobs)==CONFIG['expected_glaciers'], 'Candidate count changed')
    ref_jobs=read_csv(ROOT/'reference/glacier_worklist.csv')
    require([j['rgi_id'] for j in jobs]==[j['rgi_id'] for j in ref_jobs], 'Glacier inventory differs from historical reference')
    # Geometry versions may change decimal roundoff: compare scientific fields,
    # not serialization. A different release or membership classification fails.
    for j,r in zip(jobs,ref_jobs):
        for key in ('screen_class','membership','release_grid_row','release_grid_col','release_status'):
            require(str(j[key])==r[key], f'{j["rgi_id"]}: {key} differs')
        require(abs(float(j['polygon_overlap_besham_km2'])-float(r['polygon_overlap_besham_km2']))<=3e-6, f'{j["rgi_id"]}: polygon overlap differs')
    reference_overlap=read_csv(ROOT/'reference/glacier_overlap.csv')
    indexed=lambda rows: {(r['rgi_id'],int(r['grid_row']),int(r['grid_col'])): float(r['glacier_polygon_area_in_cell_m2']) for r in rows}
    actual, original=indexed(overlaps), indexed(reference_overlap)
    require(actual.keys()==original.keys(),'Polygon-to-cell membership differs')
    require(all(abs(actual[k]-original[k])<=max(1.0,1e-7*original[k]) for k in original),'Polygon-to-cell overlap differs')
    BASIN.parent.mkdir(parents=True,exist_ok=True)
    require(not BASIN.exists(), f'Existing basin output: {BASIN}; preserve or use separate project')
    stage=Path(tempfile.mkdtemp(prefix='.basin-',dir=BASIN.parent))
    try:
        write_csv(stage/'cells.csv',cells)
        write_csv(stage/'routing.csv',[{**v,'cell_lat':float(lat[k[0]]),'cell_lon':float(lon[k[1]])} for k,v in sorted(routes.items())])
        write_csv(stage/'candidates.csv',candidates);write_csv(stage/'overlap.csv',overlaps)
        write_csv(stage/'release_screen.csv',screen);write_csv(stage/'worklist.csv',jobs)
        with Dataset(stage/'mask.nc','w',format='NETCDF4_CLASSIC') as ds:
            ds.createDimension('lat',360);ds.createDimension('lon',720)
            ds.createVariable('lat','f8',('lat',))[:]=lat
            ds.createVariable('lon','f8',('lon',))[:]=lon
            ds.createVariable('mask','i1',('lat','lon'))[:]=mask
        atomic_json(stage/'report.json',dict(status='PROVISIONAL_BASIN_PREPARED',outlet=list(OUTLET),cells=len(cells),area_km2=float(area[mask].sum()/1e6),glaciers=len(jobs),classification=counts,rgi_sources=sources,routing_source_sha256=sha(ldd_path),area_source_sha256=sha(area_path),outputs={p.name:sha(p) for p in stage.iterdir()},limitations=['Gauge equivalence and glacier release locations are provisional.']))
        stage.rename(BASIN)
    finally:
        if stage.exists():shutil.rmtree(stage)
    print('Basin prepared:',BASIN,flush=True)


def read_jobs(root=ROOT):
    report=read_json(BASIN/'report.json')
    for name in ('worklist.csv','overlap.csv','cells.csv','routing.csv','mask.nc'):
        require(sha(BASIN/name)==report['outputs'][name], f'Basin output changed: {name}')
    jobs=read_csv(BASIN/'worklist.csv')
    require(len(jobs)==CONFIG['expected_glaciers'] and len({j['rgi_id'] for j in jobs})==len(jobs),'Worklist count or IDs invalid')
    return jobs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch-rgi',action='store_true')
    args=parser.parse_args()
    if BASIN.exists():
        read_jobs(); print('Verified existing basin:',BASIN); return
    prepare(args.fetch_rgi)

if __name__=='__main__':main()
