# Predict from an fSCAD fit

Predict from an fSCAD fit

## Usage

``` r
# S3 method for class 'fscad'
predict(object, newdata, ...)

# S3 method for class 'fscad'
fitted(object, ...)

# S3 method for class 'fscad'
residuals(object, ...)
```

## Arguments

- object:

  An object of class `fscad`.

- newdata:

  Optional finite numeric matrix whose columns correspond exactly to the
  training `argvals`. If omitted, training fitted values are returned.

- ...:

  Unused.

## Value

Numeric vector of predicted responses.
