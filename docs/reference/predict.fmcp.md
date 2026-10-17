# Predict from an fMCP fit

Predict from an fMCP fit

## Usage

``` r
# S3 method for class 'fmcp'
predict(object, newdata, ...)

# S3 method for class 'fmcp'
fitted(object, ...)

# S3 method for class 'fmcp'
residuals(object, ...)
```

## Arguments

- object:

  An object of class `fmcp`.

- newdata:

  Optional finite numeric matrix whose columns correspond exactly to the
  training `argvals`. If omitted, training fitted values are returned.

- ...:

  Unused.

## Value

Numeric vector of predicted responses.
