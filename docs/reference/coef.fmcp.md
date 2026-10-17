# Model coefficients for an fMCP fit

Model coefficients for an fMCP fit

## Usage

``` r
# S3 method for class 'fmcp'
coef(object, type = c("basis", "function"), ...)
```

## Arguments

- object:

  An object of class `fmcp`.

- type:

  Return B-spline coefficients (`"basis"`) or coefficient-function
  values on the training grid (`"function"`).

- ...:

  Unused.

## Value

A named numeric vector for `type = "basis"`, or a data frame with
`argvals` and `beta` for `type = "function"`.
