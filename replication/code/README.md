# Simulation and data-analysis code

Use R 4.1 or later, Python 3.9 or later, and a C++17 compiler configured for R. Start in
this directory after cloning the repository.

## Setup

```bash
cd replication/code
bash run.sh verify
mkdir -p library
Rscript --vanilla scripts/install_dependencies.R
R CMD INSTALL --library=library study/PFLR-master
R CMD INSTALL --library=library study/FMCP-package
R CMD INSTALL --library=library ../..
Rscript --vanilla scripts/preflight.R library
```

Package versions are listed in `software_versions.csv`. The setup script
installs dependencies into the local `library` directory.

## Prepare a run

```bash
cd replication/code
bash run.sh prepare run_01
bash run.sh acquire run_01
```

`acquire` retrieves the raw data from Zenodo and CRAN and checks their hashes.
Sources and licenses are listed in [DATA_SOURCES.md](DATA_SOURCES.md).
Each run has its own directory under `work`.

## Real-data inputs

To reproduce the inputs from the raw data:

```bash
cd replication/code
bash run.sh prepare-data run_01
mkdir -p logs
bash run.sh baselines run_01 > logs/baselines_run_01.log 2>&1
bash run.sh design run_01
```

`prepare-data` runs the original Truck/Tecator and NeoSpectra preparation
scripts. The earlier Truck pipeline also prepares Tecator; its data are
provided by `caret`. `baselines` runs the original full-data tuning,
repeated-CV and nested-CV scripts to produce the folds and comparison
predictions. This is a long numerical step, not a smoke test. Seeds, tuning
rules and preprocessing remain those in the included study scripts.
`design` then checks both study configurations without fitting their models.

If the six original inputs are already available locally, use their existing
study directory instead of `prepare-data` and `baselines`:

```bash
cd replication/code
bash run.sh import-inputs run_01 /path/to/existing/study
bash run.sh design run_01
```

Supply an absolute path to the directory containing `real_data`. The import
checks all six files against `inputs.csv` and refuses to replace existing
files. These checksums identify the original inputs; newly prepared RDS files
also contain local paths and preparation times, so their file hashes differ.
No data or prediction files are bundled in this repository.

## Simulations

```bash
cd replication/code
mkdir -p logs
bash run.sh simulation run_01 > logs/simulation_run_01.log 2>&1
bash run.sh analyze-simulation run_01
```

The simulation has 100 replications in each of four cases. The reporting step
creates figures and tables, including the equal-validation comparison.

## Real-data analyses

```bash
cd replication/code
mkdir -p logs
bash run.sh real run_01 > logs/real_run_01.log 2>&1
bash run.sh analyze-real run_01
```

Truck and NeoSpectra use the folds and Smooth/PCR predictions produced by the
baseline scripts, or the checksum-verified original inputs imported above.

## Progress and resume

```bash
cd replication/code
tail -f logs/simulation_run_01.log
```

Use `bash run.sh resume-simulation run_01` or
`bash run.sh resume-real run_01` after an interruption. Keep the same run name,
code, package versions, and tuning settings when resuming. Completed stages
are checked before they are reused. Results and reports are written under
`work/run_01`; existing output directories are protected against overwrite.
The analysis commands check output completeness, convergence, and file hashes.
Runtime depends on the computer; no current runtime estimate is available.

To resume an interrupted baseline CV stage:

```bash
cd replication/code
CODE_ROOT="$(pwd)"
export R_LIBS_USER="$CODE_ROOT/library${R_LIBS_USER:+:$R_LIBS_USER}"
export FMCP_REAL_LIBRARY="$CODE_ROOT/library"
cd work/run_01
# Run only the command for the interrupted stage.
Rscript --vanilla real_data/03_run_repeated_cv.R --resume
Rscript --vanilla real_data/11_run_neospectra_nested_cv_v3.R --resume
```

Do not rerun completed preparation or full-data tuning stages: their outputs
are protected against replacement. The full-data tuning script has no
checkpoint; if interrupted, preserve its partial directory and restart that
stage in a new run directory.
