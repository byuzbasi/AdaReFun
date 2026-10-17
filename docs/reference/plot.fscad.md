# Plot an fSCAD coefficient function

Plot an fSCAD coefficient function

## Usage

``` r
# S3 method for class 'fscad'
plot(
  x,
  shade_null = TRUE,
  null_col = "grey90",
  xlab = "Argument",
  ylab = expression(beta(t)),
  main = "fSCAD coefficient function",
  lwd = 2,
  ...
)
```

## Arguments

- x:

  An object of class `fscad`.

- shade_null:

  Logical. Shade intervals whose fitted functional norm is numerically
  zero.

- null_col:

  Color for shaded null intervals.

- xlab, ylab, main, lwd:

  Standard plotting parameters.

- ...:

  Additional arguments passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
