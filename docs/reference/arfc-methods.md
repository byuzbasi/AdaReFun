# Methods for an ARFC fit

Methods for an ARFC fit

## Usage

``` r
# S3 method for class 'arfc'
coef(object, type = c("basis", "function"), ...)

# S3 method for class 'arfc'
predict(object, newdata, ...)

# S3 method for class 'arfc'
fitted(object, ...)

# S3 method for class 'arfc'
residuals(object, ...)

# S3 method for class 'arfc'
nobs(object, ...)

# S3 method for class 'arfc'
print(x, ...)

# S3 method for class 'arfc'
summary(object, ...)

# S3 method for class 'summary.arfc'
print(x, ...)

# S3 method for class 'arfc'
plot(
  x,
  type = c("coefficient", "regions", "weights", "tuning"),
  shade_null = TRUE,
  null_col = "grey90",
  ...
)
```

## Arguments

- object, x:

  An object of class `arfc`.

- type:

  For [`coef()`](https://rdrr.io/r/stats/coef.html), return basis
  coefficients or coefficient-function values. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), choose the
  coefficient curve, regional norms, adaptive weights, or tuning path.

- ...:

  Additional graphical arguments where applicable.

- newdata:

  Optional finite predictor matrix on the training grid.

- shade_null:

  Shade finest-level regions with numerically zero fitted regional norm
  in coefficient plots.

- null_col:

  Color used for zero-region shading.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a vector or data
frame; [`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return numeric
vectors; [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns an
integer; printing and plotting methods return their object invisibly.
