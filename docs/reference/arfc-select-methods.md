# Methods for an ARFC-Select fit

Methods for an ARFC-Select fit

## Usage

``` r
# S3 method for class 'arfc_select'
coef(object, type = c("basis", "function"), ...)

# S3 method for class 'arfc_select'
predict(object, newdata, ...)

# S3 method for class 'arfc_select'
fitted(object, ...)

# S3 method for class 'arfc_select'
residuals(object, ...)

# S3 method for class 'arfc_select'
nobs(object, ...)

# S3 method for class 'arfc_select'
print(x, ...)

# S3 method for class 'arfc_select'
summary(object, ...)

# S3 method for class 'summary.arfc_select'
print(x, ...)

# S3 method for class 'arfc_select'
plot(x, type = c("validation", "coefficient", "candidates"), ...)
```

## Arguments

- object, x:

  An object of class `arfc_select`.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), show
  validation losses, the selected coefficient, or all candidate
  coefficient functions. For
  [`coef()`](https://rdrr.io/r/stats/coef.html), passed to the selected
  [`arfc()`](https://byuzbasi.github.io/AdaReFun/reference/arfc.md) fit.

- ...:

  Additional arguments.

- newdata:

  Optional predictor matrix on the training grid.

## Value

Methods delegate numerical results to the selected `arfc` fit; plotting
and printing return `x` invisibly.
