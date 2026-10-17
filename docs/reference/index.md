# Package index

## Adaptive regional models

One specification or an independent-validation choice among four.

- [`arfc()`](https://byuzbasi.github.io/AdaReFun/reference/arfc.md) :
  Fit adaptive regional folded-concave functional regression
- [`arfc_select()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_select.md)
  : Select the ARFC penalty family and resolution by validation loss
- [`arfc_control()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_control.md)
  : Control adaptive regional folded-concave fitting

## Fixed-resolution comparators

- [`fmcp()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp.md)
  [`fMCP()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp.md)
  [`FMCP()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp.md) :
  Fit functional MCP (fMCP) regression
- [`fscad()`](https://byuzbasi.github.io/AdaReFun/reference/fscad.md)
  [`fSCAD()`](https://byuzbasi.github.io/AdaReFun/reference/fscad.md) :
  Fit functional SCAD regression
- [`fmcp_control()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp_control.md)
  : Control numerical details of fMCP and fSCAD fitting

## Inspect, predict and visualize

- [`coef(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`predict(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`fitted(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`residuals(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`nobs(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`print(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`summary(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`print(`*`<summary.arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  [`plot(`*`<arfc>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-methods.md)
  : Methods for an ARFC fit
- [`coef(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`predict(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`fitted(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`residuals(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`nobs(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`print(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`summary(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`print(`*`<summary.arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  [`plot(`*`<arfc_select>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/arfc-select-methods.md)
  : Methods for an ARFC-Select fit
- [`coef(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/coef.fmcp.md)
  : Model coefficients for an fMCP fit
- [`coef(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/coef.fscad.md)
  : Model coefficients for an fSCAD fit
- [`predict(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/predict.fmcp.md)
  [`fitted(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/predict.fmcp.md)
  [`residuals(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/predict.fmcp.md)
  : Predict from an fMCP fit
- [`predict(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/predict.fscad.md)
  [`fitted(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/predict.fscad.md)
  [`residuals(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/predict.fscad.md)
  : Predict from an fSCAD fit
- [`plot(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/plot.fmcp.md)
  : Plot an fMCP coefficient function
- [`plot(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/plot.fscad.md)
  : Plot an fSCAD coefficient function
- [`print(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/print.fmcp.md)
  : Print an fMCP fit
- [`print(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/print.fscad.md)
  : Print an fSCAD fit
- [`summary(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/summary.fmcp.md)
  [`print(`*`<summary.fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/summary.fmcp.md)
  : Summarize an fMCP fit
- [`summary(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/summary.fscad.md)
  [`print(`*`<summary.fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/summary.fscad.md)
  : Summarize an fSCAD fit
- [`nobs(`*`<fmcp>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/nobs.fmcp.md)
  : Number of observations in an fMCP fit
- [`nobs(`*`<fscad>`*`)`](https://byuzbasi.github.io/AdaReFun/reference/nobs.fscad.md)
  : Number of observations in an fSCAD fit

## Package

- [`AdaReFun`](https://byuzbasi.github.io/AdaReFun/reference/AdaReFun-package.md)
  [`AdaReFun-package`](https://byuzbasi.github.io/AdaReFun/reference/AdaReFun-package.md)
  : AdaReFun: Adaptive Regional Functional Regression
