# Weak-instrument-robust confidence sets for proxy-SVAR impulse responses

Montiel Olea, Stock & Watson (2021) Anderson-Rubin-type confidence sets
for the impulse responses to a shock identified with a single external
instrument. For each variable \\i\\ and horizon \\h\\, the set collects
all \\\lambda\\ for which the (delta-method, HAC) Wald test of \\s\\
e_i'\Phi_h\Gamma - \lambda\Gamma_j = 0\\ does not reject, where \\\Gamma
= E\[u_t m_t\]\\, \\j\\ is the `target` variable and \\s\\ is `size`.
The test is a quadratic \\a\lambda^2 + b\lambda + c \le 0\\, so the set
is one of four types (reported in `case`):

- 1:

  `a > 0`, discriminant `> 0`: bounded interval `[lower, upper]`.

- 2:

  `a < 0`, discriminant `> 0`: the *complement* of an interval, i.e.
  `(-Inf, r1] U [r2, Inf)`. `lower = -Inf`, `upper = Inf` and the
  excluded interval `(r1, r2)` is in `gap`.

- 3:

  `a > 0`, discriminant `< 0`: empty set (`lower = upper = NA`); signals
  misspecification.

- 4:

  otherwise (`a < 0`, discriminant `<= 0`): the whole real line (`-Inf`,
  `Inf`); the instrument is uninformative.

- 0:

  impact response of `target`, normalised to `size`.

Unlike the bootstrap bands of
[`varm_proxy()`](https://kvasilopoulos.github.io/varm/reference/varm_proxy.md),
these have correct coverage regardless of instrument strength. The
covariance of \\(\mathrm{vec}(A), \Gamma)\\ uses the OLS influence
functions with a Bartlett (Newey-West) HAC estimator; derivatives of the
MA matrices follow Luetkepohl (2005).

## Usage

``` r
varm_proxy_robust(
  fit,
  m,
  target = 1L,
  size = 1,
  horizon = 20L,
  level = c(0.68, 0.9),
  nwlags = NULL
)
```

## Arguments

- fit:

  a [`varm()`](https://kvasilopoulos.github.io/varm/reference/varm.md)
  object.

- m:

  proxy: numeric vector, or `T x r` matrix of `r` proxies, aligned with
  the data (the first `p` rows are dropped to match the residuals).

- target:

  variable(s) (name or position) whose impact response is normalised to
  `size`; length `r`.

- size:

  size of the shock (impact response of `target`); scalar or length `r`.

- horizon:

  maximum horizon (impact is 0).

- level:

  confidence levels (in `(0, 1)`).

- nwlags:

  Newey-West truncation lag; default `floor(4 * (T/100)^(2/9))`.

## Value

a `varm_proxy_robust` object: `estimate` (`(H+1) x K` point IRF),
`lower`, `upper`, `case` (arrays `(H+1) x K x length(level)`), `gap`
(`(H+1) x K x length(level) x 2`, excluded interval in case 2), `fstat`
(see
[`proxy_fstat()`](https://kvasilopoulos.github.io/varm/reference/proxy_fstat.md)),
`level`.

## References

Montiel Olea, J. L., Stock, J. H. and Watson, M. W. (2021). Inference in
structural vector autoregressions identified with an external
instrument. *Journal of Econometrics* 225(1), 74-87.

## Examples

``` r
set.seed(1)
e <- matrix(rnorm(600), 300)
y <- matrix(0, 300, 2)
for (t in 2:300) y[t, ] <- 0.5 * y[t - 1, ] + c(e[t, 1] + 0.3 * e[t, 2], e[t, 2] + 0.5 * e[t, 1])
fit <- varm(y, p = 1)
r <- varm_proxy_robust(fit, m = e[, 1] + rnorm(300, sd = 0.5), horizon = 6)
r$lower[, , 2]
#>             y1          y2
#> 0  1.000000000  0.44467847
#> 1  0.347960448  0.10377790
#> 2  0.111524006 -0.01139702
#> 3  0.030156304 -0.03282190
#> 4  0.004596911 -0.02779278
#> 5 -0.001986987 -0.01862251
#> 6 -0.002744429 -0.01125363
```
