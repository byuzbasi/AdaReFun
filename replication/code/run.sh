#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export ARFC_FMCP_LIBRARY="$ROOT/library"
export ARFC_PFLR_LIBRARY="$ROOT/library"
export FMCP_REAL_LIBRARY="$ROOT/library"
export R_LIBS_USER="$ROOT/library${R_LIBS_USER:+:$R_LIBS_USER}"
mode="${1:-help}"
case "$mode" in
  verify)
    python3 "$ROOT/scripts/verify_manifest.py"
    ;;
  prepare)
    python3 "$ROOT/scripts/prepare_run.py" "${2:?Provide a new run name}"
    ;;
  acquire|import-inputs|prepare-data|baselines|design|simulation|real|resume-simulation|resume-real|analyze-simulation|analyze-real)
    run="${2:?Provide the prepared run name}"
    [[ "$run" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "Invalid run name"; exit 2; }
    cd "$ROOT/work/$run"
    if [ "$mode" = acquire ]; then
      Rscript --vanilla "$ROOT/scripts/acquire_raw.R"
      exit 0
    fi
    if [ "$mode" = import-inputs ]; then
      python3 "$ROOT/scripts/import_inputs.py" "${3:?Provide the existing study directory}"
      exit 0
    fi
    Rscript --vanilla "$ROOT/scripts/preflight.R" "$ROOT/library"
    case "$mode" in
      prepare-data)
        Rscript --vanilla real_data/01_prepare_real_data.R
        Rscript --vanilla real_data/10_prepare_neospectra_v3.R
        ;;
      baselines)
        Rscript --vanilla real_data/02_fit_full_models.R
        Rscript --vanilla real_data/03_run_repeated_cv.R
        Rscript --vanilla real_data/11_run_neospectra_nested_cv_v3.R
        ;;
      design)
        Rscript --vanilla simulation/45_run_arfc_select_v1.R --validate-design
        Rscript --vanilla real_data/26_run_arfc_real_data_v1.R --validate-design
        ;;
      simulation|resume-simulation)
        extra=(); [ "$mode" = resume-simulation ] && extra=(--resume)
        for stage in adaptive legacy; do
          if [ "$stage" = adaptive ]; then
            runner=simulation/45_run_arfc_select_v1.R
            completed=simulation/results/arfc_select_v1/production_run
          else
            runner=simulation/49_run_arfc_legacy_benchmark_v1.R
            completed=simulation/results/arfc_legacy_benchmark_v1/production_run_v1
          fi
          if [ "$mode" = resume-simulation ] && [ -d "$completed" ]; then
            Rscript --vanilla "$ROOT/scripts/verify_stage.R" "$completed"
          else
            Rscript --vanilla "$runner" "${extra[@]}"
          fi
        done
        ;;
      real|resume-real)
        extra=(); [ "$mode" = resume-real ] && extra=(--resume)
        Rscript --vanilla real_data/26_run_arfc_real_data_v1.R "${extra[@]}"
        ;;
      analyze-simulation)
        Rscript --vanilla simulation/46_analyze_arfc_select_v1.R
        Rscript --vanilla simulation/47_write_arfc_select_v1.R --report-version=v2
        Rscript --vanilla simulation/48_validate_arfc_select_v1.R --report-version=v2
        Rscript --vanilla simulation/50_analyze_arfc_legacy_benchmark_v1.R
        Rscript --vanilla simulation/51_write_arfc_legacy_benchmark_v1.R --report-version=v2
        Rscript --vanilla simulation/52_validate_arfc_legacy_benchmark_v1.R --report-version=v2
        Rscript --vanilla simulation/59_refine_arfc_legacy_figures_v8.R
        Rscript --vanilla "$ROOT/scripts/recompute_ise.R" config/generated/v18
        Rscript --vanilla analysis/equal_validation_v1_20261016/analyze.R
        Rscript --vanilla analysis/equal_validation_v1_20261016/verify.R
        ;;
      analyze-real)
        Rscript --vanilla real_data/27_analyze_arfc_real_data_v1.R
        Rscript --vanilla real_data/28_write_arfc_real_data_v1.R --report-version=v4
        Rscript --vanilla real_data/29_validate_arfc_real_data_v1.R --report-version=v4
        ;;
    esac
    ;;
  *)
    echo "Usage: bash run.sh verify | prepare NAME | acquire NAME | design NAME"
    echo "Inputs: prepare-data NAME | baselines NAME | import-inputs NAME EXISTING_STUDY"
    echo "Run: simulation NAME | real NAME"
    echo "Resume: resume-simulation NAME | resume-real NAME"
    echo "Report and check: analyze-simulation NAME | analyze-real NAME"
    ;;
esac
