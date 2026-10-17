# Model coefficients for an fSCAD fit

Model coefficients for an fSCAD fit

## Usage

``` r
# S3 method for class 'fscad'
coef(object, type = c("basis", "function"), ...)
```

## Arguments

- object:

  An object of class `fscad`.

- type:

  Return B-spline coefficients (`"basis"`) or coefficient-function
  values on the training grid (`"function"`).

- ...:

  Unused.

## Value

A named numeric vector for `type = "basis"`, or a data frame with
`argvals` and `beta` for `type = "function"`.
