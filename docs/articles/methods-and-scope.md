# Methods, selection and inferential scope

## Which function should I use?

| Function | Role | Selection |
|:---|:---|:---|
| arfc() | One ASR/AMR and MCP/SCAD specification | Training BIC on supplied grids |
| arfc_select() | Choose one of the four adaptive specifications | Supplied independent validation |
| fMCP() | Fixed-resolution nonadaptive MCP comparator | Training BIC |
| fSCAD() | Fixed-resolution nonadaptive SCAD comparator | Training BIC |

ARFC is the **adaptive regional folded-concave family**. ASR uses one
finest partition; AMR combines multiple partitions of the same domain,
sharing one spline coefficient vector. ARFC-Select chooses a penalty
family and resolution using response prediction error, not truth-aware
support metrics.

## What is normalized?

Let the functional domain be `[0,T]`, let `I` be an interval of positive
length, and let the coefficient curve be `beta`. A regional magnitude is
its RMS value on `I`, and its penalty multiplier is the interval’s
fraction of the domain. Level weights sum to one. Consequently adding a
resolution level does not multiply the accumulated regional curvature
merely by adding regions.

Pilot weights adjust regional shrinkage. The public implementation clips
them for numerical stability. MCP shape 2.7 and SCAD shape 3.7 have
equal maximum scalar concavity; this does not make their estimates
identical or establish the superiority of either penalty.

## Theory and software are related, but not interchangeable

The paper proves a mass identity, a matrix weak-convexity envelope,
idealized majorization descent, conditional regional-oracle results, and
a finite-library validation-risk bound. Oracle recovery requires
signal/score separation and design conditions. These results do not
certify every numerical initialization, post-selection confidence
intervals, or support-optimal validation selection.

The software uses a regional-norm floor, coefficient cutoff and
active-set refit. These conventions are exposed and retained in
diagnostics. Failed paths are not silently discarded from diagnostic
tables. Inspect the chosen fit’s convergence before reporting its
predictions or structural zeros.

## Sources and contribution

SCAD: Fan and Li (2001),
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273).
MCP: Zhang (2010),
[doi:10.1214/09-AOS729](https://doi.org/10.1214/09-AOS729). Functional
local sparsity: Lin et al. (2017),
[doi:10.1080/10618600.2016.1195273](https://doi.org/10.1080/10618600.2016.1195273).

The companion study concerns scale-normalized adaptive regional geometry
and joint penalty-resolution selection, not the invention of functional
MCP. Use `citation("AdaReFun")` to cite the software version used.
