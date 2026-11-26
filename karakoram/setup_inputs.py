"""Stage model inputs or capture their provenance; never delete original files.

A fresh rebuild uses raw ERA5 and RGI inputs. The optional replay shortcut
links verified prepared forcing to test CWatM without repeating glacier work.
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import shutil
import tempfile
import re
from urllib.parse import urlsplit, urlunsplit
from pathlib import Path
from .common import ROOT, CONFIG, SPATIAL, BASIN, CLIMATE, MAPS, require, sha, atomic_json, read_json, read_csv

STATIC_ROOT=Path('data/cwatm/spatial/30min')
STATIC_DIRS=('cwatm_input_30min','official_CWatM_Earth_30min')


def git_value(path,*args):
    return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()


def static_hashes(root):
    """Hash the scientific spatial inputs, excluding Git's internal storage."""
    output={}
    for folder in STATIC_DIRS:
        base=root/STATIC_ROOT/folder
        require(base.is_dir(),f'Missing CWatM spatial data: {base}')
        for current,dirs,files in os.walk(base):
            dirs[:]=sorted(d for d in dirs if d!='.git')
            for name in sorted(files):
                path=Path(current)/name
                if path.is_file() and path.suffix.lower() in ('.nc','.nc4','.map','.tif','.tiff','.xml','.csv','.asc'):
                    output[str(path.relative_to(root))]=sha(path)
    require(output,'No static maps found')
    return output


def link(source,target):
    """Link existing inputs only; never replace a pre-existing target."""
    require(source.exists(),f'Missing source: {source}')
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.is_symlink() and target.resolve()==source.resolve():return
    require(not target.exists() and not target.is_symlink(),f'Target already exists: {target}')
    target.symlink_to(source.resolve(),target_is_directory=source.is_dir())


def stage(source,reuse_forcing=False,reuse_raw=False):
    """Use an existing installation as a source without importing old scripts."""
    require(source.resolve()!=ROOT,'Source and destination must be different projects')
    for folder in STATIC_DIRS:
        link(source/STATIC_ROOT/folder,ROOT/STATIC_ROOT/folder)
    repo=ROOT/'cwatm/source'
    if not repo.exists():
        repo.parent.mkdir(parents=True,exist_ok=True)
        # Separate checkout; the running model does not use the old repository.
        subprocess.run(['git','clone','--local','--no-hardlinks',str(source/'cwatm/source'),str(repo)],check=True)
    require(git_value(repo,'rev-parse','HEAD')==CONFIG['cwatm_revision'],'CWatM checkout differs from reference')
    if reuse_raw:
        link(source/'data/climate/era5/raw',ROOT/'data/climate/era5/raw')
        # RGI source polygons only, not the derived inventory tables.
        link(source/'.runtime-home/OGGM/rgi/RGIV62',ROOT/'data/raw/rgi')
    if reuse_forcing:
        link(source/'data/climate/era5/processed/daily_2005_2019',CLIMATE)
        link(source/'data/cwatm/coupling/era5_besham_2005_2019',MAPS)
        BASIN.mkdir(parents=True,exist_ok=True)
        link(source/'data/basins/upper_indus/provisional_besham_30min_mask.nc',BASIN/'mask.nc')
    atomic_json(ROOT/'config/local_inputs.json',dict(mode='prepared-forcing replay' if reuse_forcing else 'raw-input rebuild',source=str(source),reused_raw=reuse_raw,static_sha256=static_hashes(ROOT)))
    print('Inputs staged. Prepared forcing replay is not an end-to-end rebuild.' if reuse_forcing else 'Raw inputs staged; run basin and climate preparation next.')


def optional_git_value(path,*args):
    try:
        return git_value(path,*args)
    except subprocess.CalledProcessError:
        return None


def public_remote(path):
    """Record public HTTPS locations without credentials, queries or fragments."""
    value=optional_git_value(path,'remote','get-url','origin')
    if not value: return None
    parsed=urlsplit(value)
    if parsed.scheme=='https' and parsed.hostname:
        return urlunsplit(('https',parsed.hostname,parsed.path,'',''))
    return None


def capture(source,output):
    """Record exact local map bytes and source revisions for later fresh installation."""
    repo=source/'cwatm/source';maps=source/STATIC_ROOT/'official_CWatM_Earth_30min'
    archive=source/STATIC_ROOT/'cwatm_input_30min.zip'
    result=dict(cwatm_revision=git_value(repo,'rev-parse','HEAD'),cwatm_remote=public_remote(repo),static_revision=optional_git_value(maps,'rev-parse','HEAD'),static_remote=public_remote(maps),original_spatial_archive_sha256=sha(archive) if archive.is_file() else None,static_sha256=static_hashes(source))
    atomic_json(output,result);print('Captured:',output)



def owned_target(path):
    """Reject links in destination paths: a fresh project owns its inputs."""
    path=Path(path)
    require(path.is_relative_to(ROOT),f'Destination outside project: {path}')
    for part in (path,*path.parents):
        if part==ROOT:break
        require(not part.is_symlink(),f'Destination is a symbolic link: {part}; use a fresh project directory')
    return path


def copy_file(source,target):
    """Copy immutable bytes, allowing safe restart without replacing differences."""
    source=Path(source);target=owned_target(Path(target))
    require(source.is_file(),f'Missing source: {source}')
    expected=sha(source)
    if target.exists():
        require(target.is_file() and sha(target)==expected,f'Destination differs: {target}')
        return expected
    target.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='.'+target.name+'-',dir=target.parent)
    os.close(fd)
    try:
        # GNU cp uses copy-on-write when available; otherwise makes a full copy.
        subprocess.run(['cp','--reflink=auto','--preserve=mode,timestamps','--',str(source),name],check=True)
        require(sha(name)==expected,f'Copy verification failed: {target}')
        os.replace(name,target)
    finally:Path(name).unlink(missing_ok=True)
    return expected


def copy_tree(source,target):
    source=Path(source);target=owned_target(Path(target))
    require(source.is_dir(),f'Missing input directory: {source}')
    shutil.copytree(source,target,symlinks=False,dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('.git','__pycache__'),copy_function=copy_file)


def capture_environment_locks(source):
    """Record package builds and pip-only versions without moving environments."""
    conda=source/'miniforge3/bin/conda'
    require(conda.is_file(),f'Original Conda executable missing: {conda}')
    lockdir=ROOT/'config/locks';lockdir.mkdir(parents=True,exist_ok=True)
    for name in ('oggm_env','cwatm_env'):
        prefix=source/'miniforge3/envs'/name
        command=[str(conda),'list','--prefix',str(prefix)]
        explicit=subprocess.check_output(command+['--explicit','--md5'],text=True)
        require('@EXPLICIT' in explicit,'Conda explicit export missing marker')
        # Do not silently create locks containing local package paths or private URLs.
        for line in explicit.splitlines():
            if not line or line.startswith('#') or line=='@EXPLICIT':continue
            parsed=urlsplit(line)
            require(parsed.scheme=='https' and parsed.hostname and not parsed.username and not parsed.password and not parsed.query and '/t/' not in parsed.path,
                    'Non-public or local package URL in explicit lock; inspect original environment manually')
        entries=json.loads(subprocess.check_output(command+['--json'],text=True))
        pip=[]
        for entry in entries:
            if entry.get('channel')=='pypi' or entry.get('build_string')=='pypi_0' or entry.get('build')=='pypi_0':
                package,version=entry['name'],entry['version']
                require(re.fullmatch(r'[A-Za-z0-9_.-]+',package) and re.fullmatch(r'[A-Za-z0-9_.!+\-]+',version),'Nonstandard pip requirement; inspect manually')
                pip.append(f'{package}=={version}')
        (lockdir/f'{name}.conda.txt').write_text(explicit)
        (lockdir/f'{name}.pip.txt').write_text(''.join(v+'\n' for v in sorted(pip)))
        versions=subprocess.check_output([str(prefix/'bin/python'),'-c',
            "import importlib.metadata as m,json;print(json.dumps({d.metadata['Name']:d.version for d in m.distributions()},sort_keys=True))"],text=True)
        atomic_json(lockdir/f'{name}.versions.json',json.loads(versions))
        print('Captured package specifications:',name,flush=True)


def migrate(source,copy_oggm=True):
    """Own a fresh copy of original inputs; do not stage derived model products."""
    source=source.resolve()
    require(source!=ROOT and not ROOT.is_relative_to(source),'Fresh project must be separate from original')
    # Reject the earlier prepared-input replay directory, even if its links are broken.
    for target in (ROOT/STATIC_ROOT,ROOT/'data/climate/era5/raw',ROOT/'data/raw/rgi',ROOT/'data/processed/basin',CLIMATE,MAPS,ROOT/'results'):
        owned_target(target)
    require(not any(p.exists() or p.is_symlink() for p in (BASIN,CLIMATE,MAPS,ROOT/'results')),
            'Derived outputs already exist; extract the release into a new empty project directory')
    require(not git_value(source/'cwatm/source','status','--porcelain','--untracked-files=no'),'Original CWatM tracked source has changes')
    capture(source,ROOT/'config/source_manifest.json')
    manifest=read_json(ROOT/'config/source_manifest.json')
    require(manifest['cwatm_revision']==CONFIG['cwatm_revision'],'Original CWatM revision differs from reference')
    capture_environment_locks(source)
    print('Copying original CWatM maps',flush=True)
    for folder in STATIC_DIRS:copy_tree(source/STATIC_ROOT/folder,ROOT/STATIC_ROOT/folder)
    original_archive=source/STATIC_ROOT/'cwatm_input_30min.zip'
    if original_archive.is_file():copy_file(original_archive,ROOT/STATIC_ROOT/original_archive.name)
    require(static_hashes(ROOT)==manifest['static_sha256'],'Static input hashes differ')
    print('Copying WAPDA original observations and metadata',flush=True)
    copy_tree(source/'data/observations/wapda/original',ROOT/'data/observations/wapda/original')
    station=ROOT/'data/observations/wapda/original/WAPDA/Cleaned_Data/SWHP_discharge/Indus_at_Besham_Qila_daily_discharge.csv'
    metadata=ROOT/'data/observations/wapda/original/WAPDA/Cleaned_Data/metadata/stations.csv'
    require(sha(station)==CONFIG['observation_source_sha256'],'WAPDA discharge source differs')
    require(sha(metadata)==CONFIG['observation_metadata_source_sha256'],'WAPDA station metadata source differs')
    print('Copying exact reference ERA5 raw acquisition',flush=True)
    raw=Path('data/climate/era5/raw');original=read_json(ROOT/'reference/original_era5_acquisition_manifest.json')
    plan=raw/'ecmwf_arco'/original['plan'];complete=source/plan/'COMPLETE.json'
    require(read_json(complete)==original,'Original raw acquisition manifest differs from reference')
    for entry in original['reused_cds']:
        rel=Path(entry['file'])
        require(not rel.is_absolute() and rel.is_relative_to(raw) and '..' not in rel.parts,'Invalid ERA5 source path')
        require(copy_file(source/rel,ROOT/rel)==entry['sha256'],f'CDS source differs: {rel}')
        request=(source/rel).with_suffix('.request.json')
        if request.is_file():copy_file(request,(ROOT/rel).with_suffix('.request.json'))
    for entry in original['tiles']:
        rel=Path(entry['file']);require(rel.name==str(rel),'Invalid tile filename')
        require(copy_file(source/plan/rel,ROOT/plan/rel)==entry['sha256'],f'ERA5 tile differs: {rel}')
    for name in ('acquisition_plan.json','COMPLETE.json'):
        src=source/plan/name
        if src.is_file():copy_file(src,ROOT/plan/name)
    print('Copying RGI 6.2 source polygons for regions 13, 14 and 15',flush=True)
    rgi=source/'.runtime-home/OGGM/rgi/RGIV62'
    for region in CONFIG['rgi_regions']:
        dirs=sorted(p for p in rgi.glob(region+'_rgi62_*') if p.is_dir())
        require(len(dirs)==1,f'RGI region missing or ambiguous: {region}')
        copy_tree(dirs[0],ROOT/'data/raw/rgi'/dirs[0].name)
    archives={}
    if copy_oggm:
        print('Copying pinned ERA5 L3 source TAR archives (may be large)',flush=True)
        cache=Path('.runtime-home/OGGM/download_cache')/('cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L3-L5_files/2025.6/elev_bands/ERA5/per_glacier_spinup/RGI62/b_080/L3')
        jobs=read_json(ROOT/'config/project.json')['reference_omissions']
        from .common import read_csv
        prefixes={j['rgi_id'][:11] for j in read_csv(ROOT/'reference/glacier_worklist.csv') if j['rgi_id'] not in jobs}
        for prefix in sorted(prefixes):
            rel=cache/prefix[:8]/(prefix+'.tar')
            require((source/rel).is_file(),f'Missing original archive {rel}; rerun with --download-glacier-archives to reacquire instead')
            archives[str(rel)]=copy_file(source/rel,ROOT/rel)
    repo=owned_target(ROOT/'cwatm/source')
    if not repo.exists():
        repo.parent.mkdir(parents=True,exist_ok=True)
        subprocess.run(['git','clone','--local','--no-hardlinks',str(source/'cwatm/source'),str(repo)],check=True)
        subprocess.run(['git','-C',str(repo),'checkout','--detach',CONFIG['cwatm_revision']],check=True)
        if manifest.get('cwatm_remote'):
            subprocess.run(['git','-C',str(repo),'remote','set-url','origin',manifest['cwatm_remote']],check=True)
        else:subprocess.run(['git','-C',str(repo),'remote','remove','origin'],check=True)
    require(git_value(repo,'rev-parse','HEAD')==CONFIG['cwatm_revision'],'CWatM checkout differs')
    require(not git_value(repo,'status','--porcelain','--untracked-files=no'),'CWatM checkout modified')
    require(not (repo/'.git/objects/info/alternates').exists(),'CWatM clone references external Git objects')
    atomic_json(ROOT/'config/raw_migration_report.json',dict(mode='independent raw-input rebuild',static_maps_verified=True,era5_plan=original['plan'],oggm_archive_sha256=archives,station_sha256=sha(station),metadata_sha256=sha(metadata),cwatm_revision=CONFIG['cwatm_revision']))
    print('Original inputs copied and verified. No old environments or derived outputs migrated.',flush=True)


def install_environments():
    """Install new environments into this project's own Miniforge prefix."""
    mf=owned_target(ROOT/'miniforge3');conda=mf/'bin/conda'
    require(conda.is_file(),'Install a fresh Miniforge at miniforge3/ first')
    require(not os.environ.get('KARAKORAM_MINIFORGE') or Path(os.environ['KARAKORAM_MINIFORGE']).resolve()==mf,'Unset KARAKORAM_MINIFORGE before installing')
    for name in ('oggm_env','cwatm_env'):
        prefix=owned_target(mf/'envs'/name)
        lockdir=ROOT/'config/locks';explicit=lockdir/f'{name}.conda.txt';pip=lockdir/f'{name}.pip.txt'
        require(explicit.is_file() and pip.is_file(),f'Capture package locks before installing {name}')
        env={**os.environ,'CONDA_PKGS_DIRS':str(mf/'pkgs'),'CONDA_ENVS_PATH':str(mf/'envs')}
        env.pop('PYTHONPATH',None)
        if not prefix.exists():
            subprocess.run([str(conda),'create','--yes','--prefix',str(prefix),'--file',str(explicit)],check=True,env=env)
        else:print('Checking existing new environment:',prefix,flush=True)
        def installed_versions():
            return json.loads(subprocess.check_output([str(prefix/'bin/python'),'-c',
                "import importlib.metadata as m,json;print(json.dumps({d.metadata['Name']:d.version for d in m.distributions()},sort_keys=True))"],text=True,env=env))
        expected=read_json(lockdir/f'{name}.versions.json')
        actual=installed_versions()
        if actual!=expected and pip.read_text().strip():
            subprocess.run([str(prefix/'bin/python'),'-m','pip','install','--no-deps','--no-cache-dir','-r',str(pip)],check=True,env=env)
            actual=installed_versions()
        differences={k:[v,actual.get(k)] for k,v in expected.items() if actual.get(k)!=v}
        extra=sorted(set(actual)-set(expected))
        require(not differences and not extra,f'Python package versions differ in {name}: {differences}; extra distributions: {extra}')
        print('Original Python package versions verified:',name,flush=True)



def verify_observations():
    """Check the two original user-supplied files before acquiring public inputs."""
    base=ROOT/'data/observations/wapda/original/WAPDA/Cleaned_Data'
    files={'SWHP_discharge/Indus_at_Besham_Qila_daily_discharge.csv':'observation_source_sha256',
           'metadata/stations.csv':'observation_metadata_source_sha256'}
    for rel,key in files.items():
        path=owned_target(base/rel)
        require(path.is_file(),f'Place original WAPDA file at {path}')
        require(sha(path)==CONFIG[key],f'Original WAPDA checksum differs: {path}')
    print('Original WAPDA station discharge and metadata verified.',flush=True)


def bootstrap():
    """Install a project-local Miniforge and locked model environments if needed."""
    import platform
    require(platform.system()=='Linux' and platform.machine()=='x86_64',
            'The captured environment locks target Linux x86_64')
    mf=owned_target(ROOT/'miniforge3')
    if not (mf/'bin/conda').is_file():
        require(not mf.exists(),f'Incomplete Miniforge installation: {mf}; inspect before retrying')
        release='26.9.4-2'
        name=f'Miniforge3-{release}-Linux-x86_64.sh'
        installer=ROOT/'data/raw/software'/name
        verified_download(f'https://github.com/conda-forge/miniforge/releases/download/{release}/{name}',installer,
                          'bbaf16b3897da24c4be50315934ba74848f617b4495982841d6a0a60e0edbbdf')
        subprocess.run(['bash',str(installer),'-b','-p',str(mf)],check=True)
    install_environments()


def static_subset(directory, folder, manifest):
    """Verify the exact reference file set for one static input directory."""
    prefix=str(STATIC_ROOT/folder)+'/'
    expected={k[len(prefix):]:v for k,v in manifest['static_sha256'].items() if k.startswith(prefix)}
    require(expected, f'No checksums for {folder}')
    require(directory.is_dir(), f'Missing map directory: {directory}')
    actual={str(p.relative_to(directory)):sha(p) for p in directory.rglob('*')
            if p.is_file() and '.git' not in p.relative_to(directory).parts
            and p.suffix.lower() in ('.nc','.nc4','.map','.tif','.tiff','.xml','.csv','.asc')}
    missing=sorted(set(expected)-set(actual)); extra=sorted(set(actual)-set(expected))
    changed=sorted(k for k in expected.keys() & actual.keys() if expected[k]!=actual[k])
    require(not (missing or extra or changed), f'Static map mismatch in {folder}: missing={missing}; extra={extra}; changed={changed}')
    return len(expected)


def verified_download(url, destination, expected):
    """Use local bytes first; promote a new download only after SHA-256 verification."""
    from urllib.request import urlopen
    destination=owned_target(destination)
    if destination.exists():
        require(sha(destination)==expected, f'Existing archive checksum differs: {destination}')
        return
    destination.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='.'+destination.name+'-',dir=destination.parent)
    try:
        with os.fdopen(fd,'wb') as output, urlopen(url,timeout=120) as response:
            shutil.copyfileobj(response,output,1024*1024)
        require(sha(name)==expected,f'Download checksum differs: {url}; original source must be identified, not substituted')
        os.replace(name,destination)
    finally:Path(name).unlink(missing_ok=True)


def extract_static_archive(archive, destination, manifest):
    """Extract into private staging, reject traversal, links and differing map bytes."""
    import zipfile, stat
    destination=owned_target(destination)
    require(not destination.exists(),f'Existing destination: {destination}')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.spatial-',dir=destination.parent) as tmp:
        stage=Path(tmp)
        with zipfile.ZipFile(archive) as z:
            for member in z.infolist():
                path=Path(member.filename)
                require(not path.is_absolute() and '..' not in path.parts and not stat.S_ISLNK(member.external_attr>>16),f'Unsafe ZIP member: {member.filename}')
            z.extractall(stage)
        candidate=stage/destination.name
        static_subset(candidate,destination.name,manifest)
        candidate.rename(destination)


def pinned_checkout(remote, revision, destination, verify=None):
    """Fetch a pinned public Git revision, including its Git LFS objects."""
    destination=owned_target(destination)
    if destination.exists():
        require(git_value(destination,'rev-parse','HEAD')==revision,f'Wrong source revision: {destination}')
        require(not git_value(destination,'status','--porcelain','--untracked-files=no'),f'Tracked source changes: {destination}')
        if verify:verify(destination)
        return
    destination.parent.mkdir(parents=True,exist_ok=True)
    require(shutil.which('git-lfs'),'Install Git LFS before acquiring CWatM source maps (Arch: sudo pacman -S git-lfs)')
    with tempfile.TemporaryDirectory(prefix='.checkout-',dir=destination.parent) as tmp:
        stage=Path(tmp)/'repository'
        env={**os.environ,'GIT_LFS_SKIP_SMUDGE':'1','GIT_TERMINAL_PROMPT':'0'}
        subprocess.run(['git','clone','--no-checkout',remote,str(stage)],check=True,env=env)
        subprocess.run(['git','-C',str(stage),'lfs','install','--local'],check=True)
        subprocess.run(['git','-C',str(stage),'checkout','--detach',revision],check=True,env=env)
        subprocess.run(['git','-C',str(stage),'lfs','pull','origin'],check=True)
        require(git_value(stage,'rev-parse','HEAD')==revision,'Pinned checkout failed')
        require(not git_value(stage,'status','--porcelain','--untracked-files=no'),'Downloaded tracked source modified')
        if verify:verify(stage)
        stage.rename(destination)


def acquire_static():
    """Install exact model and map sources without accessing the retired project."""
    manifest=read_json(ROOT/'config/source_manifest.json')
    require(manifest['cwatm_revision']==CONFIG['cwatm_revision'],'Source manifest revision differs')
    base=owned_target(ROOT/STATIC_ROOT)
    archive=base/'cwatm_input_30min.zip'
    legacy=base/STATIC_DIRS[0]
    if not legacy.exists():
        print('Acquire original CWatM spatial archive if missing',flush=True)
        verified_download('https://zenodo.org/records/3528097/files/cwatm_input_30min.zip?download=1',archive,manifest['original_spatial_archive_sha256'])
        extract_static_archive(archive,legacy,manifest)
    count=static_subset(legacy,STATIC_DIRS[0],manifest)
    official=base/STATIC_DIRS[1]
    if not official.exists():
        print('Acquire pinned official CWatM spatial repository',flush=True)
        pinned_checkout(manifest['static_remote'],manifest['static_revision'],official,lambda p:static_subset(p,STATIC_DIRS[1],manifest))
    count+=static_subset(official,STATIC_DIRS[1],manifest)
    pinned_checkout(manifest['cwatm_remote'],manifest['cwatm_revision'],ROOT/'cwatm/source')
    require(static_hashes(ROOT)==manifest['static_sha256'],'Static input verification failed')
    print(f'Local CWatM source and {count} spatial files verified.',flush=True)


def configure_oggm(internet):
    """Write only project-owned OGGM paths; never depend on a user's old home."""
    home=ROOT/'.runtime-home'; base=home/'OGGM'
    path=owned_target(home/'.oggm_config')
    home.mkdir(parents=True,exist_ok=True)
    values={'dl_cache_dir':base/'download_cache','dl_cache_readonly':False,
            'tmp_dir':base/'tmp','rgi_dir':base/'rgi','test_dir':base/'tests',
            'has_internet':bool(internet)}
    fd,name=tempfile.mkstemp(prefix='.oggm_config-',dir=home)
    try:
        with os.fdopen(fd,'w') as f:
            for key,value in values.items():f.write(f'{key} = {value}\n')
        os.replace(name,path)
    finally:Path(name).unlink(missing_ok=True)


def missing_rgi_regions():
    """Locate complete RGI shapefile sets in either supported local source."""
    missing=[]
    for region in CONFIG['rgi_regions']:
        found=False
        for base in (ROOT/'data/raw/rgi',ROOT/'.runtime-home/OGGM/rgi'):
            paths=sorted(base.rglob(f'{region}_rgi62_*.shp'))
            require(len(paths)<=1,f'Ambiguous RGI source: {region}')
            if paths:
                require(all(paths[0].with_suffix(ext).is_file() for ext in ('.shp','.dbf','.shx','.prj')),
                        f'Incomplete local RGI source: {paths[0]}')
                found=True;break
        if not found:missing.append(region)
    return missing


def glacier_archives(reproduce_omissions):
    """Use regenerated worklist, not reference output, to select archive groups."""
    worklist=read_csv(ROOT/'data/processed/basin/worklist.csv')
    require(len(worklist)==CONFIG['expected_glaciers'],'Unexpected glacier worklist count')
    omissions=set(CONFIG['reference_omissions']) if reproduce_omissions else set()
    url=CONFIG['era5_archive'].rstrip('/')+'/RGI62/b_080/L3/'
    prefixes=sorted({j['rgi_id'][:11] for j in worklist if j['rgi_id'] not in omissions})
    result={}
    for group in prefixes:
        require(re.fullmatch(r'RGI60-(13|14)\.\d{2}',group),'Invalid archive group')
        remote=url+group[:8]+'/'+group+'.tar'
        u=urlsplit(remote)
        local=ROOT/'.runtime-home/OGGM/download_cache'/u.netloc/u.path.lstrip('/')
        result[remote]=owned_target(local)
    return result


def auxiliary_present():
    """Known sample-data caches for the pinned OGGM/Salem installation."""
    home=ROOT/'.runtime-home'
    salem=home/'.salem_cache/salem-sample-data-454bf696324000d198f574a1bf5bc56e3e489051'
    oggm=home/'.oggm/oggm-sample-data-8af40f89620c6bd72f3485a777a018dcacb99d94'
    return salem.is_dir() and any(salem.iterdir()) and oggm.is_dir() and any(oggm.iterdir())


def prepare_oggm(stage,reproduce_omissions=False,allow_download=False):
    """Acquire missing source groups once; computation always uses offline mode.

    Does not touch glacier calculations, the worklist or numerical provenance.
    Run with the isolated OGGM launcher. Network errors here stop the workflow
    before the glacier loop, rather than marking individual glaciers blocked.
    """
    regions=missing_rgi_regions() if stage=='basin' else []
    archives=glacier_archives(reproduce_omissions) if stage=='glaciers' else {}
    # Missing or zero-length archives need attention. Existing nonempty TARs
    # are checked by OGGM when extracted; corrupt local files are not silently replaced.
    missing={url:path for url,path in archives.items() if not path.is_file() or not path.stat().st_size}
    auxiliary=stage=='glaciers' and not auxiliary_present()
    needs_network=bool(regions or missing or auxiliary)
    configure_oggm(False)
    print(f'OGGM input check: {len(regions)} missing RGI regions, '
          f'{len(missing)}/{len(archives)} missing glacier archives, '
          f'auxiliary initialization needed: {auxiliary}',flush=True)
    if not needs_network:
        print('Local inputs available; OGGM offline mode enabled.',flush=True)
        return
    require(allow_download,'Missing OGGM source inputs; explicitly allow acquisition with --download')
    configure_oggm(True)
    try:
        import oggm
        from oggm import cfg,utils
        require(oggm.__version__==CONFIG['oggm_version'],'Unexpected OGGM version')
        # Full initialization acquires the auxiliary data for a cold installation.
        cfg.initialize(logging_level='WARNING')
        for region in regions:
            print('Acquire RGI region',region,flush=True)
            utils.get_rgi_region_file(region,version='62')
        for url,path in missing.items():
            print('Acquire glacier archive',path.name,flush=True)
            actual=utils.file_downloader(url)
            require(actual and Path(actual).resolve()==path.resolve() and path.stat().st_size>0,
                    f'Archive acquisition failed: {url}')
        require(not (missing_rgi_regions() if stage=='basin' else []),'RGI acquisition incomplete')
        require(all(p.is_file() and p.stat().st_size>0 for p in archives.values()),'Glacier archive acquisition incomplete')
    finally:
        configure_oggm(False)
    print('Acquisition complete; OGGM offline mode enabled.',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
    a=sub.add_parser('stage');a.add_argument('source',type=Path);a.add_argument('--reuse-forcing',action='store_true');a.add_argument('--reuse-raw',action='store_true')
    a=sub.add_parser('capture');a.add_argument('source',type=Path);a.add_argument('--output',type=Path,default=ROOT/'config/source_manifest.json')
    a=sub.add_parser('migrate');a.add_argument('source',type=Path);a.add_argument('--download-glacier-archives',action='store_true',help='Do not copy cached OGGM archives; download them later')
    sub.add_parser('install-environments')
    sub.add_parser('acquire-static')
    sub.add_parser('bootstrap')
    sub.add_parser('verify-observations')
    a=sub.add_parser('prepare-oggm');a.add_argument('stage',choices=['basin','glaciers']);a.add_argument('--download',action='store_true');a.add_argument('--reproduce-omissions',action='store_true')
    a=p.parse_args()
    if a.mode=='stage':stage(a.source.resolve(),a.reuse_forcing,a.reuse_raw)
    elif a.mode=='capture':capture(a.source.resolve(),a.output)
    elif a.mode=='migrate':migrate(a.source,not a.download_glacier_archives)
    elif a.mode=='prepare-oggm':prepare_oggm(a.stage,a.reproduce_omissions,a.download)
    elif a.mode=='acquire-static':acquire_static()
    elif a.mode=='bootstrap':bootstrap()
    elif a.mode=='verify-observations':verify_observations()
    else:install_environments()

if __name__=='__main__':main()
