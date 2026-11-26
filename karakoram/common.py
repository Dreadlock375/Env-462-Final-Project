"""Project layout, shared checks and atomic output helpers."""
from pathlib import Path
import csv
import hashlib
import json
import os
import tempfile
from datetime import date, timedelta

ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / 'config/project.json').read_text())
SPATIAL = ROOT / 'data/cwatm/spatial/30min/cwatm_input_30min'
BASIN = ROOT / 'data/processed/basin'
CLIMATE = ROOT / 'data/climate/era5/processed/daily_2005_2019'
GLACIERS = ROOT / 'results/glaciers/historical'
MAPS = ROOT / 'data/processed/glacier_forcing'
HYDROLOGY = ROOT / 'results/hydrology'
OBS = ROOT / 'data/processed/observations'
START = date.fromisoformat(CONFIG['start'])
END = date.fromisoformat(CONFIG['end'])
DAYS = [START + timedelta(days=i) for i in range((END-START).days+1)]
OUTLET = tuple(CONFIG['outlet_row_col'])


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.'+path.name+'-', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(obj, f, indent=2, allow_nan=False)
            f.write('\n')
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    require(bool(rows), f'No rows to write: {path}')
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
