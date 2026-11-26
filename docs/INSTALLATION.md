# Installation and acquisition

## Independent rebuild

The release contains the exact captured source manifest, Linux-64 Conda explicit package locks (including package MD5 values), pip-only version lists and installed Python distribution versions. No previous installation is needed. Pip versions are locked, but wheel hashes and system libraries are not locked; cross-machine bitwise reproducibility is not guaranteed. Git and Git LFS must be available on the host. The pinned Miniforge 26.9.4-2 Linux x86_64 installer is downloaded from the official release and verified against its published SHA-256. An existing project-owned Miniforge is reused. Environment creation uses the explicit package lists, and distribution versions are verified. An incomplete installation is not overwritten automatically.

Provide your original WAPDA files in `data/observations/wapda/original/WAPDA/Cleaned_Data/` at the paths listed in the root README. If required, place your personal CDS credentials in `.runtime-home/.cdsapirc` and accept the ERA5 dataset licence. `CDSAPI_RC` may explicitly select another credentials file. The workflow automatically sets HOME and cache paths to this project; no old environment override is supported for the independent workflow. Host libraries and network access remain external dependencies.

```bash
bash workflow.sh inputs bootstrap
bash workflow.sh inputs acquire-static
bash workflow.sh test
bash workflow.sh rebuild-source
```

Or run `rebuild-source` directly after supplying the observations and credentials. The three setup commands above are optional checkpoints.

## Static inputs

The legacy spatial archive candidate is the public CWatM release at https://zenodo.org/records/3528097. Before extracting, verify its SHA-256 against `config/source_manifest.json`. The archive must contain the `cwatm_input_30min/` directory. The second spatial map set is obtained at the recorded revision of `https://github.com/iiasa/CWatM-Earth-30min.git` with Git LFS, and CWatM at the recorded revision of `https://github.com/iiasa/CWatM.git`. Both static map sets must match all 141 captured checksums exactly. Existing correct files are used without contacting the server. Existing invalid files cause an error; the workflow never silently replaces them. The actual remote archive bytes and Git checkouts have not been downloaded/tested in the development environment. If an upstream source differs, preserve the original verified bytes and resolve the source rather than accepting changed maps.

## Climate and glacier sources

A complete local ERA5 acquisition is reused. When missing, the existing ECMWF ARCO acquisition module obtains the required hours and checksums. The exact original mixed CDS/ARCO acquisition can be preserved for strict reference reproduction; a new all-ARCO acquisition is not guaranteed byte identical. OGGM acquisition checks the regenerated worklist for all required archive groups. Missing RGI polygons, archives, or required auxiliary data are acquired first, then OGGM computation uses local caches with `has_internet=False`. Connection failures should stop acquisition, not silently produce glacier omissions. Existing completed scientific outputs are validated by their respective stages. Acquisition requires availability of the public upstream repositories and authorized CDS access. An existing incomplete nonempty archive must be investigated rather than overwritten.

## Current-project upgrade

Apply the release only after active processes stop. It does not change `project.json`, the numerical glacier driver, physics, basin, coupling, climate or hydrology modules relative to the v2/local-first implementation. The source manifest and lock files are captured from the user's actual configuration. Existing datasets are not included or overwritten. `inputs acquire-static` is a useful local verification command for the current migration. Run `test` and retain the rebuilt validation and regression reports. `rebuild-reference` requires all 5,113 daily calibrated outputs to match the reference within rtol=1e-8, atol=1e-5 m³/s. `rebuild-source` records regression differences for newly acquired data; it is not permission to claim an exact match. Both modes repeat the original calibration search and retain the existing calibration-selection gate.
