# recoverR

**Calibrated layer attribution and spatial recoverability for bulk multi-omics.**

`recoverR` answers two questions that existing multi-omics tooling does not.

---

## 1. Does adding an omics layer actually improve out-of-sample prediction?

Not "does it explain variance in a fitted model" — that is a different, easier,
and much less useful question.

* `variancePartition` decomposes variance across **named covariates** within a fitted model.
* `MOFA2` decomposes variance across **latent factors** within a fitted model.

Both are descriptive decompositions of a model that has already seen all the data.
Neither reports an **out-of-sample** attribution, and neither reports when an
added layer is actively *harmful*.

`recoverR` computes the **exact Shapley (LMG) decomposition** of a
cross-validated performance functional, plus a **dimension-matched permutation
null**, so "this layer adds nothing" comes back as a p-value and a confidence
interval rather than a bare R².

## 2. Can a bulk profile carry a spatially-defined phenotype at all?

Reliability layers now exist for **spatial deconvolution output** — they add a
risk map to RCTD/cell2location/Tangram predictions and ask *is this spot-level
spatial call trustworthy?* They all operate spatial → spatial.

None asks the inverse question, which is the one that decides whether a bulk
biomarker programme is viable: **given a bulk profile, is the spatial archetype
recoverable — and if not, can the method say so?**

`recoverR` provides split and Mondrian conformal prediction with an explicit
**abstention** rule. Abstention is the deliverable, not a shortcoming: a method
that returns a confident archetype for a profile that cannot support one is the
mechanism by which bulk biomarkers fail in clinical trials.

---

## Core algorithms

Throughout, **larger is better** for every functional. The **full** specification — assumptions,
validity arguments, complexity, and the failure mode that forced each non-negotiable design choice — is
in **[`ALGORITHMS.md`](ALGORITHMS.md)**; this is the summary.

### 1. The value functional — out-of-sample performance of a set of views $S$

$$
\mathcal{V}(S) \;=\; \mathcal{P}\Big( y,\ \{\, \hat f^{(S)}_{-k}(x_i) \,\}_{i \in F_k,\ k = 1..m} \Big)
$$

$\hat f^{(S)}_{-k}$ is trained on every fold but $k$, and $\mathcal{P}$ is the **pooled** functional
(pooling rather than averaging per-fold scores, because per-fold AUC is undefined when a fold holds few
positives). **One fold vector is generated once and reused for every subset $S$** — with independent
random splits, split noise at the $n \in [20, 300]$ typical of these cohorts dominates the effect being
measured, and every subset is then trained on identical folds. Within each training fold, all
standardisation and hyper-parameter selection happens **inside that fold only** (nested selection).

$\mathcal{P}$ is one of:

$$
R^2_{\mathrm{oos}} \;=\; 1 - \frac{\sum_{i \in \text{test}} (y_i - \hat y_i)^2}{\sum_{i \in \text{test}} (y_i - \bar y_{\text{test}})^2}
\qquad\text{(Gaussian)}
$$

$$
\mathrm{AUC} \;=\; \frac{\sum_{i:\,y_i=1} r_i \;-\; n_1(n_1+1)/2}{n_1 n_0}
\qquad\text{(binomial; $r_i$ = midrank of $\hat p_i$)}
$$

$$
C \;=\; \frac{\#\{(i,j) : \eta_i > \eta_j\} + \tfrac{1}{2}\#\{(i,j) : \eta_i = \eta_j\}}{\#\{(i,j)\ \text{comparable}\}}
\qquad\text{(survival; $\eta$ = risk score)}
$$

The AUC is computed from **midranks** rather than through `pROC::roc.test()`, which errors at
$\mathrm{AUC} = 1$ — precisely where the strongest biomarkers sit — and returns `NA`, that then
propagates silently through any downstream false-discovery adjustment.

### 2. Layer attribution — the exact Shapley (LMG) decomposition

With the views as players and $\mathcal{V}$ as the payoff:

$$
\phi_k \;=\; \sum_{S \subseteq V \setminus \{k\}} \frac{|S|!\,(K - |S| - 1)!}{K!}\,\Big[\mathcal{V}(S \cup \{k\}) - \mathcal{V}(S)\Big]
$$

This is the unique attribution satisfying **efficiency** ($\sum_k \phi_k = \mathcal{V}(V) - \mathcal{V}(\varnothing)$),
**symmetry**, **dummy** and **additivity** — so a redundant layer is not credited for correlation with a
useful one, and the answer cannot depend on the order the views happened to be listed. Cost is
$2^K \times C_{\mathrm{cv}}$; above $K = 6$ the package **stops** rather than approximating.

Two derived readouts:

$$
U_k \;=\; \mathcal{V}(V) - \mathcal{V}(V \setminus \{k\}) \quad\text{(the "unique" contribution)},
\qquad
\text{forward: } k^{(t)} = \arg\max_{k \notin S_{t-1}}\Big[\mathcal{V}(S_{t-1} \cup \{k\}) - \mathcal{V}(S_{t-1})\Big]
$$

* $U_k \approx 0$ with $\phi_k > 0$: the layer carries information already present elsewhere — **redundant**, not worthless.
* $U_k < 0$: removing the layer **improves** held-out performance — the layer **actively harms**. No variance decomposition can see this; it is the strongest form of the package's central finding.

### 3. The dimension-matched permutation null

$$
X_k^{(\pi)} = P_\pi X_k,
\qquad
\Delta_k = \mathcal{V}(V) - \mathcal{V}(V \setminus \{k\}),
\qquad
\hat p = \frac{1 + \#\{\,b : \Delta_k^{(b)} \ge \Delta_k \,\}}{1 + B_{\mathrm{finite}}}
$$

Permuting view $k$'s **rows** destroys its sample-to-outcome association while preserving its internal
covariance $\Sigma_k$, its marginal distributions, its sparsity — and **exactly $p_k$ features**, so the
null is dimension-matched **by construction** and needs no dimensionality correction. The hypothesis is
precise:

> $H_0$: given the other views, view $k$ carries no sample-specific information about $y$.

The $+1$ in numerator and denominator is the finite-$B$ correction: it makes the test valid rather than
merely asymptotically valid, and prevents a reported $\hat p = 0$, which is never a legitimate result of a
finite permutation test. The null **mean** is reported too, because a layer whose $\Delta$ falls *below*
it is contributing less than a randomly re-paired version of itself.

### 4. Recoverability — split conformal prediction with abstention

Nonconformity scores $s_i$ on a calibration set of size $n_c$, the finite-sample quantile, and the
prediction set:

$$
\hat q = s_{(\lceil (n_c+1)(1-\alpha) \rceil)}, \qquad s_{(n_c+1)} := +\infty,
\qquad
C(x) = \{\, y : s(x,y) \le \hat q \,\}
$$

which yields **exact finite-$n$ coverage** under exchangeability alone (no distributional form):

$$
\mathbb{P}\big(Y_{n+1} \in C(X_{n+1})\big) \;\ge\; 1 - \alpha
$$

The $\lceil (n_c+1)(1-\alpha)\rceil$ correction is what delivers finiteness — the plain empirical quantile
under-covers at the $n \approx 20$–$50$ that immuno-oncology cohorts actually have. That guarantee is
**marginal**; **Mondrian** (class-conditional) calibration supplies the per-group version, calibrating a
separate quantile per archetype $g$ and admitting each candidate label against **its own** quantile:

$$
\hat q_g = s^{g}_{(\lceil (n_g+1)(1-\alpha) \rceil)},
\qquad
\mathbb{P}\big(Y \in C(X) \ \big|\ Y = g\big) \;\ge\; 1-\alpha \quad \forall g
$$

**Abstention is the deliverable**, not a shortcoming:

$$
\mathrm{abstain}(x) = \mathbb{1}\big[\,|C(x)| > 1\,\big]
$$

a *singleton* set is a committed call; a larger set is the method declining to commit. And the actionable
output is the **confusability** matrix (**high = a bulk profile cannot separate the pair**):

$$
\Sigma_{ij} = \frac{1}{n_c}\sum_{t=1}^{n_c} \mathbb{1}\big[\{i,j\} \subseteq C_t\big],
\qquad i \ne j
$$

so the paper-style statement is *which clinical distinction a bulk assay will fail to deliver* rather than
one aggregate accuracy that hides it.

### 5. Sampling adequacy — how much tissue, and how many genes?

Regions of diameter $d$ are pseudo-bulked from spot counts ($x_R = \sum_{s \in R} c_s$, or the mean), the
region label being the majority spot label. Performance against diameter is then fitted by nonlinear least
squares with a **Hill** curve:

$$
P(d) \;=\; P_\infty \, \frac{d^{\,h}}{K^{\,h} + d^{\,h}}
$$

$P_\infty$ is the plateau (asymptotically recoverable signal), $K$ the half-saturation diameter and $h$ the
cooperativity. A Hill form is used rather than a single exponential because it accommodates a genuine
**threshold** (large $h$) as well as a gradual rise ($h \approx 1$); a one-parameter exponential forces a
shape and therefore biases the answer. Solving $P(d^\star) = t\,P_\infty$:

$$
d^\star \;=\; K \left(\frac{t}{1-t}\right)^{1/h}
$$

$d^\star$ is the **minimum region diameter** at which the endpoint is recoverable at fraction $t$ (default
$0.95$) of its asymptote — a measured, data-driven **design rule** for how much tissue and how many genes a
study needs. Where the fitted $d^\star$ lies beyond the sampled range the package sets
`extrapolated = TRUE` and you report *"the required diameter exceeds what we measured"* rather than a
number.

### The full specification

| § of [`ALGORITHMS.md`](ALGORITHMS.md) | what it adds beyond the above |
|---|---|
| 0 | notation |
| 1 | the three functionals, with their ranges and degeneracies |
| 2 | the estimator, **nested selection**, the shared-fold argument, the `type.measure` trap |
| 3 | Shapley vs a variance share (the four axioms), the unique/shared split, exact enumeration |
| 4 | what the null preserves and destroys, the finite-$B$ correction, `null_q95` |
| 5 | leave-one-out and forward selection read together |
| 6 | the LAC/APS scores, weighted conformal under covariate shift, reliability & ECE, the full abstention/selective-risk curve, and where recoverability is deliberately **not** invented |
| 7 | pseudo-bulk construction, the Hill fit's bound checks, why the scorer is deliberately simple |
| 8–10 | cost per step, **what these algorithms do not claim**, and the underlying references |

## Install

**Requirements:** R **≥ 4.1.0**. The seven `Imports` — `ggplot2`, `matrixStats`, `Matrix`, `MASS`,
`stats`, `utils`, `grDevices` — install automatically.

### From GitHub (recommended)

```r
# remotes — the lightest option
install.packages("remotes")
remotes::install_github("LingzhangMeng/recoverR", build_vignettes = TRUE)

# or pak — faster, resolves dependencies in parallel
install.packages("pak")
pak::pak("LingzhangMeng/recoverR")

# or devtools
install.packages("devtools")
devtools::install_github("LingzhangMeng/recoverR")
```

### From a clone or a downloaded tarball

```bash
git clone https://github.com/LingzhangMeng/recoverR.git
cd recoverR
R CMD INSTALL .                 # install the checkout
R CMD build .                   # …or build recoverR_0.1.3.tar.gz
R CMD check recoverR_0.1.3.tar.gz
```

```r
# in R, from the checkout directory:
pak::local_install(".")          # or  devtools::install(".")
# …or from a built tarball (a PATH — see the note below):
install.packages("/path/to/recoverR_0.1.3.tar.gz", repos = NULL, type = "source")
```

> **⚠️ A common mistake.** `install.packages("recoverR", repos = NULL, type = "source")` — a bare
> **package name** — does not install anything: with `repos = NULL` the first argument must be a **file
> path** (a tarball) or you should use `R CMD INSTALL <dir>` / `pak::local_install("<dir>")` for a
> directory. Use the GitHub route above unless you are working from a checkout.

### If GitHub is slow or blocked (China)

```bash
git clone https://ghproxy.net/https://github.com/LingzhangMeng/recoverR.git
cd recoverR && R CMD INSTALL .
```
Any working mirror of the repository URL can be substituted; the package itself has no download-time
dependency on GitHub.

### Optional — the packages that unlock specific features

Everything below is in `Suggests`. `recoverR` installs and runs without them, and tells you when a feature
needs one (`need_pkg()`).

| what you want | install |
|---|---|
| the cross-validated learner used by default | `install.packages("glmnet")` |
| `Surv` endpoints (Harrell's $C$) | `install.packages("survival")` |
| the AUC helpers | `install.packages("pROC")` |
| extra plot layers (plate, curves) | `install.packages(c("patchwork", "cowplot", "ggrepel", "ggdist", "ggridges", "viridis", "RColorBrewer", "scales"))` |
| `read_visium()` (10x `.h5` counts) | `BiocManager::install("rhdf5")` |
| heatmap-style figures | `BiocManager::install(c("ComplexHeatmap", "circlize"))` |
| SVG/PNG output, font handling | `install.packages(c("svglite", "ragg", "systemfonts"))` |
| the MOFA+/SNF wrappers | `BiocManager::install(c("MOFA2", "SNFtool"))` |
| `data.table` / `jsonlite` helpers | `install.packages(c("data.table", "jsonlite"))` |
| running the test suite | `install.packages("testthat")` |

`BiocManager` itself, if you do not have it: `install.packages("BiocManager")`.

### Verify

```r
library(recoverR)
packageVersion("recoverR")     # → 0.1.3
vignette("recoverR")           # the package vignette
?attrib_views                  # the top-level entry point
```

## Tutorial

**New here? Read [`TUTORIAL.md`](TUTORIAL.md).** It is a complete runnable walkthrough — every module
with its **real printed output** and figures, how to read each figure, the guardrails, what is and is
not validated, and the traps that quietly produce plausible wrong numbers. Its figures are regenerated
by `Rscript tutorial/make_tutorial.R` (deterministic; seeds are listed there).

![Attribution plate: Shapley shares, leave-one-view-out loss, forward selection](man/figures/tut-02-plate.jpeg)

## Quick start

```r
library(recoverR)

# X_list: named list of view matrices, samples x features, rows aligned
# y: a two-level factor, a numeric vector, or a Surv object
res <- attrib_views(X_list, y, family = "binomial", n_perm = 199)
res
#> <recoverR_attribution> 4 views | binomial endpoint
#>
#> Shapley decomposition of out-of-sample performance:
#>    view   shapley  V_alone  unique_contribution
#>     RNA    0.4021   0.5502               0.3180
#>     ...
#>
#> Dimension-matched permutation null:
#>   Methylation    delta = +0.0041   p = 0.6200  (B = 199)
#>   ...

plot_shapley(res)                     # hero panel
plot_attribution_plate(res)           # Shapley + LOO + forward
save_plot(plot_shapley(res), "figs/shapley")   # writes .pdf AND .jpeg
```

### Calibrated recoverability

```r
fit <- refit_recoverability(bulk_X, archetype, alpha = 0.10, mondrian = TRUE)
pred <- predict_recoverability(fit, new_bulk_X)

table(pred$abstain)                    # who the model refuses to call
archetype_separability(fit)            # which pairs bulk cannot separate
plot_separability(fit)                 # the actionable figure
```

### Sampling adequacy

```r
sa <- sampling_adequacy(spot_counts, coords, labels,
                        diameters = c(0.5, 1, 2, 3, 4, 5, 6))
sa
#> Adequacy thresholds:
#>   panel_genes  plateau  half_saturation_mm  adequacy_mm
#>           200    0.842                1.85         3.42
```

---

## The guardrails

Every function in `R/guardrails.R` exists because a real analysis silently
produced a wrong answer **without raising an error**. They are exported because
the failure modes are general:

| Wrapper | Failure it prevents |
|---|---|
| `cv_glmnet_safe()` | `cv.glmnet` silently downgrades `type.measure` from `"auc"` to `"deviance"` when folds hold < 10 observations, so `which.max(cv$cvm)` maximises **error** |
| `auc_safe()` | `pROC::roc.test()` errors at AUC = 1 — exactly where the strongest markers sit — and returns `NA`, which propagates through BH adjustment |
| `safe_cbind()` | Per-layer PCs are all named `PC1..PCn`, so `cbind` + `lm` returns `NA` for the combined model |
| `assert_no_na()` | `if (x) next` over an `NA`-containing vector fails with "missing value where TRUE/FALSE needed" |
| `detect_cna_levels()` | `paste0("C", NA)` yields the literal `"CNA"`, adding a phantom factor level that breaks `scale_fill_manual` at draw time |

## Figures

`save_plot()` writes **both** a vector PDF and a JPEG, and then **verifies** the
PDF is genuinely vector. That check matters: ggplot2 ≥ 3.5 rasterises continuous
colourbar legends by default, embedding a one-pixel bitmap strip in the PDF, and
`raster = FALSE` is deprecated and silently ignored in those versions. Use
`rr_guide_colourbar()` (or a discrete fill) and the figure stays editable.

```r
save_plot(p, "figs/panel")                      # pdf + jpeg, vector-verified
save_pub(p, "figs/panel1", width_mm = 183)      # svg + pdf + 600 dpi tiff
```

Plot functions: `plot_shapley()`, `plot_leave_one_out()`,
`plot_forward_selection()`, `plot_attribution_plate()`, `plot_ablation_bars()`,
`plot_variance_partition()`, `plot_reliability()`, `plot_coverage_risk()`,
`plot_separability()`, `plot_recoverability()`, `plot_sampling_curve()`.

---

## What this does **not** claim

* Shapley attribution is **associative, not causal**.
* Conformal coverage is **marginal** unless `mondrian = TRUE`, and degrades under
  distribution shift.
* **Abstention is not accuracy.** A model that abstains on hard cases is honest,
  not better.
* Nothing here substitutes for **external validation**.

Every exported function documents its own mathematics and behaviour (`?attrib_views`,
`?refit_recoverability`, `?sampling_adequacy`), including the
validity arguments, the complexity analysis, and the optimisations.

## Citation

If you use `recoverR`, cite the package and the accompanying analysis. See
`inst/CITATION`.

## License

MIT © Lingzhang Meng
