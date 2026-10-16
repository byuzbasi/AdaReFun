# Contributing

Report issues with a minimal reproducible example, `sessionInfo()`, the
explicit tuning grids, integration rule, and relevant convergence diagnostics.
Do not upload personal, confidential, or restricted data.

Numerical changes require tests against the established reference kernels and
clear documentation of effects on scientific results. Preserve public APIs and
retain failed candidate diagnostics. Run `R CMD check` before proposing changes.

Documentation and plotting improvements should not alter fitting defaults.
Reproducing the article uses its frozen study engines, not a silent substitution
of a newer package version. The authors welcome discussion before substantial
methodological changes.
