# varm plot theme and palette

`theme_varm()` is the single theme used by every
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
method in varm: black-and-white panel, dashed major grid, no strip
background, legend at the bottom. `varm_palette()` returns colours from
the colour-blind safe Okabe-Ito palette.

## Usage

``` r
theme_varm(base_size = 11)

varm_palette(n = 8L)
```

## Arguments

- base_size:

  base font size, passed to
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

- n:

  number of colours (at most 8; more are recycled with a warning).

## Value

`theme_varm()`: a
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html);
`varm_palette()`: a character vector of hex colours.

## See also

Other plotting:
[`autoplot.varm_draws()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_draws.md),
[`autoplot.varm_forecast()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_forecast.md),
[`autoplot.varm_hd()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_hd.md),
[`autoplot.varm_irf()`](https://kvasilopoulos.github.io/varm/reference/autoplot.varm_irf.md)

## Examples

``` r
varm_palette(3)
#> [1] "#0072B2" "#D55E00" "#009E73"
ggplot2::ggplot(yraw, ggplot2::aes(infl, rate)) +
  ggplot2::geom_point(colour = varm_palette(1)) +
  theme_varm()
```
