# Data sources, access and redistribution

## NeoSpectra

Sanderman, Partida, Safanelli, Mitu, Ge, Murad and Shepherd (2024),
*Near-infrared (NIR) soil spectral library using the NeoSpectra Handheld NIR
Analyzer by Si-Ware*, version 1.2, https://doi.org/10.5281/zenodo.13122321.
The record's license is Creative Commons Attribution 4.0 (CC BY 4.0), verified
from the Zenodo record on 16 October 2026.

Original file: `1976_NSlibrary_withmetadata.csv`.
MD5: `8cb55e45af554a69bf2b818c9ff31492`.
Acquisition script uses the immutable record's file endpoint. The raw file is
not bundled. The processed RDS is a derivative: 1,974 retained sample IDs,
601 wavelengths, within-sample averaging of repeat spectra, and response
`log1p(SOC percent)` with source field `oc_usda.c729_w.pct`. The processed
object retains the source manifest and preparation metadata. The derivative
is generated locally, not bundled in this repository. It remains subject to
the source's CC BY 4.0 terms, rather than the software GPL.

Description paper: Partida et al. (2025), *Data in Brief* 58, 111229,
https://doi.org/10.1016/j.dib.2024.111229.

## Truck

The distributed Truck data are obtained from PFLR 1.1.0 (GPL-2). The PFLR
documentation attributes the data to the Coordinating Research Council
E55/E59 study, Clark et al. (2007). The source package and its data documentation
are included for attribution; raw `truck.rda` is retrieved separately from the
official source archive and checked against MD5
`597a06ef330f8ff33dc38a419676940d`.

The processed derivative retains the original 108 observations, response and
90 acceleration measurements. It adds the normalized argument grid used by the
analysis, without changing raw values. Its provenance remains attached to the
PFLR distribution; retain GPL-2 notices with redistributed derivatives. No
vehicle identifier is available; no inferred identifier has been added.

## Input roles

Processed RDS objects are fitting inputs. Fold CSVs define the same held-out
samples across methods. Smooth/PCR prediction CSVs are secondary benchmark
inputs, not training data. The original preparation and comparison scripts
generate these files locally. `inputs.csv` records the six original input
checksums and their generating scripts for users importing an existing study.
The earlier Truck pipeline also prepares the public-domain Tecator data from
`caret`; Tecator is not part of the final Truck/NeoSpectra application.
Bootstrap/outer-CV result files are outputs. No raw file is edited by
acquisition, preparation, estimation, plotting or reporting commands.
