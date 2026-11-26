# Karakoram — historical OGGM–CWatM

Self-contained historical workflow for the provisional Besham Qila Upper Indus basin, 2005–2019. This release preserves the frozen historical scientific implementation and does not require any old scripts, W5E5 data, previous project, or previous environment. Future scenarios are not implemented yet.

## Start from scratch

Requirements: Linux x86_64, system Python 3, Bash, Git, Git LFS, ordinary system libraries, sufficient disk space and access to upstream data services. On Arch install Git LFS if missing (`sudo pacman -S git-lfs`). Place the original user-supplied WAPDA directory at `data/observations/wapda/original/WAPDA/`. At minimum it must contain `Cleaned_Data/SWHP_discharge/Indus_at_Besham_Qila_daily_discharge.csv` and `Cleaned_Data/metadata/stations.csv`. Supply your own licensed CDS credentials in `.runtime-home/.cdsapirc` (mode 600). Do not commit or upload this file. The additional observational archives can be kept under `data/observations/` but are not needed by this historical model.

From the project root at its permanent absolute path:

```bash
bash workflow.sh inputs verify-observations
bash workflow.sh rebuild-source
```

This checks WAPDA input hashes, installs a pinned project-local Miniforge if absent, recreates the locked environments, fetches and verifies both static map sets and the pinned CWatM revision, reconstructs the basin and inventory, acquires hourly ERA5, converts daily climate, acquires RGI/OGGM inputs as needed, computes glacier runoff, exports CWatM forcing, repeats calibration and evaluates the frozen parameters. The source-mode workflow records any reference daily-discharge difference rather than asserting exact equality. If the calibration search chooses different parameters, the existing gate stops for review. It retains the seven reference glacier omissions, fixed historical geometry and basin assumptions for comparability.

`bash workflow.sh rebuild-reference` follows the same sequence but demands exact daily-discharge regression within the existing numerical tolerances. Keep the original mixed CDS/ARCO raw inputs if exact historical reproduction is important; fresh all-ARCO acquisition may differ from them. Neither mode requires the old project.

## Individual operations

| Command (`bash workflow.sh ...`) | Operation |
|---|---|
| `inputs bootstrap` | Create or verify local software environments |
| `inputs acquire-static` | Acquire missing CWatM inputs and check all 141 map hashes |
| `inputs verify-observations` | Verify user-supplied WAPDA source files |
| `test` | Automated numerical, integrity and regression tests |
| `basin --fetch-rgi` | Reconstruct basin and glacier worklist |
| `observations` | Prepare daily WAPDA records |
| `acquire-climate` | Acquire original ERA5 hours |
| `prepare-climate` | Generate daily meteorological forcing |
| `glaciers --fetch-prepro --prune-gdirs --reproduce-omissions` | Pre-acquire glacier archives, then calculate offline |
| `export --accept-omissions` | Generate melt/rain/cover maps |
| `historical --repeat-calibration --require-reference-match` | Repeat calibration and validation with strict regression |

See `docs/INSTALLATION.md` for details and `docs/METHODOLOGY.md` for scientific assumptions. Do not move the installed Conda prefix. Do not delete the old project until the currently rebuilt results have passed the independent regression and input checks.

## Data access and publication

This source distribution does not include WAPDA observations or the raw climate,
static-map or glacier archives. Obtain the original WAPDA station CSV and
metadata independently from the data provider or project owner, subject to the
provider's permission. This historical reference configuration requires their
exact recorded checksums; other datasets are not automatically interchangeable.
The reference daily CSV contains simulated discharge only. A SHA-256 fingerprint
of the original observation sequence preserves the exact observation regression
check without redistributing individual measured discharge values. The two
observation-dependent tests skip if the original WAPDA data are absent; supply
the observations to run them. Neither a successful test suite with skips nor a
source-only checkout establishes a complete cold rebuild.

The reference rebuild with migrated raw inputs and fresh environments reproduced
all 5,113 historical daily simulated discharge values exactly. The automated
cold acquisition workflow has not yet been fully exercised from observations
alone. The model remains a historical research workflow with provisional basin
and glacier routing, fixed glacier geometry, and imperfect flood timing and
frequency. Future climate scenarios are not yet implemented.

Third-party software and datasets retain their respective licensing and access
conditions; downloading this repository does not grant redistribution rights to
the original WAPDA records or upstream inputs.
