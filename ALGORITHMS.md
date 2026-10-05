# recoverR — Mathematical Algorithms

**Version 0.1.3 · 2026-09-16, revised 2026-09-29**

This document specifies every algorithm the package implements, with the
estimator, its assumptions, its validity argument, and its computational cost.
Where a design choice was forced by a failure mode observed in practice, the
failure mode is stated — those are the reasons the choice is not negotiable.

> **0.1.3 (2026-09-29) — the archetype rule is now shared, and the imaging route can feed it.**
> The three-way call (`f_immune < desert_cut` -> desert; else `f_tumour_immune >= inflamed_cut` ->
> inflamed; else excluded) existed only inside `derive_spatial_archetypes()`, which takes *spot-level
> programme scores*. An **immunofluorescence** study measures the same three programmes directly and
> reports them as **per-ROI counts**, so it cannot call that function — and giving it a second copy of
> the rule would let the two routes drift, invisibly, until someone compared two published archetype
> calls. The rule is therefore extracted as the internal `call_archetype()`, used by both.
>
> **New export `read_roi_table()`** shapes a segmentation's per-ROI counts
> (`n_tumour`, `n_stroma`, `n_tumour_immune`, `n_stroma_immune`) into exactly the columns the Visium
> route produces (`n_spots`, `n_tumour`, `f_immune`, `f_tumour_immune`, `f_stroma_immune`, `archetype`,
> `region`, `sample_id`), so `fit_recoverability()` -> `predict_recoverability()` ->
> `reliability_curve()` -> `archetype_separability()` -> `sampling_adequacy()` run on imaging data
> **unchanged**. It fails loudly on the traps that matter: a **fraction or area supplied where a count
> belongs** (caught as "more immune cells than compartment cells"), a missing or non-numeric column, and
> negative counts. Regions with **no tumour compartment** are reported as `NA` — undefined rather than
> forced into a class — and counted in an attribute. Extra columns pass through, so a ROI's size can
> travel with it, which is what a spatial-scale question needs. Five tests pin the equivalence to the
> Visium route and the guards (`tests/testthat/test-spatial-io.R`).

> **0.1.2 (2026-09-28) — the LDA ill-conditioning criterion, corrected.** `learn_probs()`'s
> auto-switch fired only at `p >= n`. That is **not** sufficient for an estimable pooled covariance:
> a real bulk-to-spatial call runs the conformal layer with **p = 2000 panel genes and n_train = 2235**, so
> the guard stayed silent and `MASS::lda()` inverted a 2000 x 2000 covariance from 2235 samples,
> returning posteriors that **saturate at ~1.0 while the call is at chance**. Because the conformal
> layer builds its per-class quantiles *from those posteriors*, every class was then admitted for
> every sample — which is why Module 2's real-data run reported "separability ~0.99 for every pair"
> (including the diagonal), ECE 0.475 and ~100 % abstention. Measured directly on that matrix with a
> **planted, zero-noise** signal, i.e. one that is definitely present:
>
> | learner | accuracy | mean confidence | ECE | abstention |
> |---|---|---|---|---|
> | `lda` (the old default path) | **0.725** | 0.999 | 0.274 | **0.996** |
> | `glmnet` | **0.954** | 0.926 | **0.029** | **0.000** |
>
> The criterion is now a **ratio** — at least `5` training samples per feature — which subsumes the
> old test and leaves `p << n` behaviour (`p = 6`, `n = 400`) untouched. Three regression tests pin
> it in both directions, including the planted-signal recovery above
> (`tests/testthat/test-conformal.R`), which is the positive control that found it.

> **0.1.1 (2026-09-27) — the `fit_predict()` fallback defect, fixed.** When every inner
> `cv.glmnet()` errored, the no-`lambda` fallback returned the **whole** `glmnet` path as an
> `n_test × 100` matrix; `cv_perf()`'s `oof[te] <- …` then **recycled** the first `n_test` values
> into the held-out fold, scoring it with an **undeclared** model behind a single warning.
> Measured at n = 16: one fold returned 500 predictions for 5 held-out rows and moved an arm from
> 0.6000 to 0.6545. The fallback now **declares one lambda** (the most regularised point of the
> same path — the column the recycling silently used, so the change is conservative), asserts the
> prediction's length, and `cv_perf()` asserts it again at the assignment site. Two regression
> tests pin both (see `tests/testthat/test-attribution.R`). The published analyses in this
> repository were re-run to **measure** the effect: their p-values and observed deltas are
> unchanged.

---

## 0. Notation

| Symbol | Meaning |
|---|---|
| $n$ | number of samples (patients) |
| $K$ | number of omics views (layers) |
| $X_k \in \mathbb{R}^{n \times p_k}$ | view $k$; rows are samples, aligned across $k$ |
| $V = \lbrace 1,\dots,K\rbrace $ | the full view set |
| $S \subseteq V$ | a subset of views |
| $X_S = [\ X_k\ ]_{k \in S}$ | column-wise concatenation over $S$ |
| $y$ | endpoint: $\mathbb{R}^n$ (gaussian), $\lbrace 0,1\rbrace ^n$ (binomial), or $(\text{time},\text{event})^n$ (survival) |
| $\mathcal{V}(S)$ | **out-of-sample** performance of a model trained on $X_S$ |
| $\hat q$ | conformal quantile (calibration) |
| $\alpha$ | miscoverage level; coverage target $1-\alpha$ |

**The central distinction.** `variancePartition` and `MOFA2` report a
decomposition of a *fitted* model: they answer "how much of the variance in the
data does this term account for, given that the model already saw everything?"
$\mathcal{V}(S)$ above is deliberately different — it is measured on data the
model did **not** see, and it is the quantity that decides whether a biomarker
will work in a new cohort.

---

## 1. Performance functionals

All functionals are oriented so that **larger is better**.

### 1.1 Gaussian — out-of-sample $R^2$

$$
R^2_{\text{oos}} = 1 - \frac{\sum_{i \in \text{test}} (y_i - \hat y_i)^2}{\sum_{i \in \text{test}} (y_i - \bar y_{\text{test}})^2}
$$

Computed on pooled out-of-fold predictions. Note the denominator uses the *test*
mean, which makes the statistic honest but gives it a range $(-\infty, 1]$; a
value below 0 means the model is worse than predicting the test mean.

### 1.2 Binomial — AUC via the Mann–Whitney form

$$
\text{AUC} = \frac{1}{n_1 n_0}\sum_{i: y_i = 1}\ \sum_{j: y_j = 0} \Big[ \mathbb{1}(\hat p_i > \hat p_j) + \tfrac{1}{2}\mathbb{1}(\hat p_i = \hat p_j) \Big]
= \frac{\sum_{i:y_i=1} r_i - n_1(n_1+1)/2}{n_1 n_0}
$$

where $r_i$ is the midrank of $\hat p_i$. This is computed **directly from ranks**
rather than through `pROC::roc.test()`, because that function errors when
$\text{AUC}=1$ — precisely where the strongest biomarkers sit — and returns `NA`,
which then propagates silently through any downstream false-discovery
adjustment. The rank form is exact, tie-aware, and never degenerates.

### 1.3 Survival — Harrell's concordance index

For all comparable pairs $(i,j)$ where $i$ experienced the event before $j$ was
censored:

$$
C = \frac{\lvert\lbrace (i,j) : \eta_i > \eta_j\rbrace \rvert + \tfrac{1}{2}\lvert\lbrace (i,j): \eta_i = \eta_j\rbrace \rvert }{\lvert\lbrace (i,j)\ \text{comparable}\rbrace \rvert}
$$

with $\eta$ the model's risk score (the linear predictor).

---

## 2. Cross-validated $\mathcal{V}(S)$ — and why the folds are shared

### 2.1 Estimator

For a fold assignment $\mathcal{F} = \lbrace F_1,\dots,F_m\rbrace $ (stratified on the
outcome),

$$
\mathcal{V}(S) = \mathcal{P}\Big( y,\ \lbrace \ \hat f^{(S)}_{-k}(x_i) \ \rbrace _{i \in F_k,\ k=1..m} \Big)
$$

where $\hat f^{(S)}_{-k}$ is trained on $\bigcup_{j \ne k} F_j$ and $\mathcal P$ is
the pooled functional of §1. Pooling rather than averaging per-fold scores is
required because per-fold AUC with a handful of positives is undefined in many
folds.

### 2.2 Nested selection (non-negotiable)

Within each training fold, **all** standardisation and hyper-parameter selection
happens on that fold only:

1. compute $\mu, \sigma$ from $X_S^{\text{train}}$ and apply to train **and**
   test with those frozen statistics;
2. select $(\alpha, \lambda)$ by an inner cross-validation on the training fold
   alone;
3. refit on the training fold and predict the held-out fold.

Selecting features or tuning on the full dataset is the single largest source of
optimistic bias in published multi-omics benchmarks. `recoverR` structures this
so it cannot happen by accident — see `fit_predict()`.

### 2.3 Shared fold assignment (statistical and computational)

One fold vector is generated once per analysis and reused for **every** subset
$S$:

- **Statistically:** comparing $\mathcal{V}(S)$ and $\mathcal{V}(S')$ computed on
  independent random splits confounds the subset effect with split noise. At the
  $n \in [20, 300]$ typical here, that noise dominates the effect being measured.
- **Computationally:** the fold structure is computed once, and every subset is
  trained on identical folds, so no work is duplicated.

### 2.4 The `type.measure` trap

`glmnet::cv.glmnet()` **silently downgrades** `type.measure = "auc"` to
`"deviance"` whenever any fold holds fewer than 10 observations, emitting only a
warning. `cv$cvm` then contains deviance, so the common idiom
`which.max(cv$cvm)` maximises *error*. At $n = 27$ with 5 folds — a typical
immuno-oncology cohort — this downgrade is the rule, not the exception.

`recoverR` therefore requests deviance explicitly (so behaviour does not depend
on $n$) and computes AUC itself from held-out predictions. `cv_glmnet_safe()`
exists to make the old behaviour impossible to invoke unknowingly.

---

## 3. Exact Shapley (LMG) decomposition of $\mathcal{V}$

### 3.1 Definition

$$
\phi_k = \sum_{S \subseteq V \setminus \lbrace k\rbrace } \frac{|S|!\ (K-|S|-1)!}{K!}\ \Big[\mathcal{V}(S \cup \lbrace k\rbrace ) - \mathcal{V}(S)\Big]
$$

This is the **exact** Shapley value of the cooperative game
$\mathcal{G}(S) = \mathcal{V}(S)$, where the "players" are omics views and the
"payoff" is out-of-sample performance. It coincides with the LMG (Lindeman,
Merenda, Gold) decomposition of regression $R^2$, and with $\mathcal{V}$
measured out-of-sample it becomes a *predictive* attribution.

### 3.2 Why Shapley rather than a variance share

$\phi_k$ is the unique attribution satisfying:

| Axiom | Statement | Why it matters here |
|---|---|---|
| **Efficiency** | $\sum_k \phi_k = \mathcal{V}(V) - \mathcal{V}(\varnothing)$ | The shares reconcile to the whole; nothing is double-counted or lost. |
| **Symmetry** | Equally contributive views receive equal $\phi$ | The answer does not depend on the order views were listed. |
| **Dummy** | A view contributing nothing to every $S$ receives $\phi_k = 0$ | A useless layer is not credited for correlation with a useful one. |
| **Additivity** | $\phi$ of a sum of games is the sum of $\phi$ | Results compose when endpoints or sub-cohorts are combined. |

The naive alternative — $R^2$ of layer $k$ alone, or the increment when added in
one fixed order — violates at least one axiom. Incremental-in-one-order is
especially misleading: a redundant layer added *last* looks worthless, while the
same layer added *first* looks essential. With 4 views, the ordering-dependent
answers span a factor that frequently exceeds the effect size being reported.

### 3.3 Decomposition of the value into unique and shared parts

Two derived quantities are reported:

- **Unique contribution** $U_k = \mathcal{V}(V) - \mathcal{V}(V \setminus \lbrace k\rbrace )$
  — the leave-one-out loss, a lower-bound reading of "what this layer buys".
- **Shapley value** $\phi_k$ — the fair share, which credits a layer for
  performance it shares with others.

Where $U_k \approx 0$ but $\phi_k > 0$, the layer carries information already
present elsewhere (**redundant**, not worthless). Where $U_k < 0$, removing the
layer **helps**: the layer actively harms out-of-sample performance. That case is
invisible to any variance decomposition, and is the strongest form of the
package's central finding.

### 3.4 Complexity and exact enumeration

$$
\text{cost} = \underbrace{2^{K}}_{\text{subsets}} \times \underbrace{C_{\text{cv}}}_{\text{fits per subset}}
$$

The $\mathcal{V}(S)$ values are computed once over the whole subset lattice and
**reused across all $K$ players**, so the total cost is one pass over the lattice
rather than $K$ passes. Subset keys are canonicalised in **view order** (not
alphabetically) because `combn()` and the Shapley loop emit views in different
orders; a mismatched key fails as `subscript out of bounds` rather than as a
wrong number, which is why the lattice size is asserted to equal $2^K$.

For $K \le 6$ this is $64 \times C_{\text{cv}}$ fits — tractable. Above
$K = 6$ the function **stops** rather than silently approximating; permutation
sampling of the lattice would be the correct extension and is deliberately not
implemented rather than approximated badly.

---

## 4. Dimension-matched permutation null

### 4.1 Construction

To test whether view $k$ contributes anything, its **rows are permuted**:

$$
X_k^{(\pi)} = P_\pi X_k, \qquad P_\pi\ \text{a random permutation matrix}
$$

The statistic $\Delta_k = \mathcal{V}(V) - \mathcal{V}(V\setminus\lbrace k\rbrace )$ is
recomputed on $X_k^{(\pi)}$ for $\pi_1,\dots,\pi_B$.

### 4.2 What this null preserves, and what it destroys

| | Destroyed | Preserved |
|---|---|---|
| View $k$ | the sample-to-profile association with $y$ | internal covariance $\Sigma_k$, marginal distributions, sparsity |
| Other views | nothing | their full alignment with $y$ |
| Geometry | nothing | exactly $p_k$ features — the null is dimension-matched **by construction** |

This is why the null needs no dimensionality correction. Comparing against
random noise of a different width would require one, and would still not test the
right hypothesis.

The null hypothesis is therefore precisely:

> $H_0$: given the other views, view $k$ carries no sample-specific information about $y$.

### 4.3 p-value

$$
\hat p = \frac{1 + \lvert\lbrace \ b : \Delta_k^{(b)} \ge \Delta_k \ \rbrace \rvert }{1 + B_{\text{finite}}}
$$

The $+1$ in numerator and denominator is the standard finite-$B$ correction; it
makes the test valid (not merely asymptotically valid) and prevents a reported
$\hat p = 0$, which is never a legitimate result from a finite permutation test.
$B_{\text{finite}} \le B$ allows for permutations that fail to fit.

Critical values are read from the same empirical null (`null_q95`), and the null
mean is reported because a layer with $\Delta$ **below** the null mean is
contributing less than a randomly re-paired version of itself — direct evidence
of harm.

---

## 5. Ablation: leave-one-out and forward selection

**Leave-one-out.** $\text{LOO}_k = \mathcal{V}(V) - \mathcal{V}(V\setminus\lbrace k\rbrace )$.
Positive $\Rightarrow$ the layer helps; negative $\Rightarrow$ it harms.

**Forward selection.** Starting from $\varnothing$, repeatedly add the view with
the largest marginal gain:

$$
k^{(t)} = \arg\max_{k \notin S_{t-1}} \Big[\mathcal{V}(S_{t-1} \cup \lbrace k\rbrace ) - \mathcal{V}(S_{t-1})\Big]
$$

Read together with $\phi_k$, the pair separates three distinct situations that a
single $R^2$ table conflates: a view that is **useful alone**, one that is
**useful only in combination**, and one that is **actively harmful**. The forward
path also exposes the greedy-order artefact directly, since the staircase
flattening is visible rather than summarised away.

---

## 6. Split conformal prediction with abstention

### 6.1 Assumption

Calibration and test points are **exchangeable**. This is the only assumption;
no distributional form is required. It is also the assumption that fails under
cohort shift (§6.5), which is exactly the regime of a multi-cohort biomarker, so
it is checked rather than assumed.

### 6.2 Algorithm

1. Split the data into a proper training set ($n_t$) and a calibration set ($n_c$).
2. Fit $\hat p = \hat f(X_{\text{train}})$.
3. Compute nonconformity scores $s_i$ on the calibration set (§6.3).
4. Compute the finite-sample quantile

$$
\hat q = s_{(\lceil (n_c+1)(1-\alpha)\rceil)}, \qquad s_{(n_c+1)} := +\infty
$$

5. Predict the set $C(x) = \lbrace y : s(x,y) \le \hat q \rbrace $.

### 6.3 Nonconformity scores

**LAC** (least ambiguous set) — threshold on the true-class probability:

$$
s_{\text{LAC}}(x,y) = 1 - \hat p_y(x)
\quad\Longrightarrow\quad
C(x) = \lbrace \ y : \hat p_y(x) \ge 1 - \hat q \ \rbrace 
$$

**APS** (adaptive prediction sets) — cumulative mass of labels at least as likely
as the true one:

$$
s_{\text{APS}}(x,y) = \sum_{j\ :\ \hat p_j(x) \ \ge\ \hat p_y(x)} \hat p_j(x)
\quad\Longrightarrow\quad
C(x) = \text{smallest top-mass prefix with cumulative} \ \ge \hat q
$$

APS is the default because it **adapts set size to difficulty**: easy samples get
singleton sets and hard samples are permitted a larger set. For this package's
purpose that is the whole point — the size of $C(x)$ *is* the recoverability
signal. LAC tends to produce either too-small or uniformly large sets and thus
carries less information about individual difficulty.

### 6.4 Coverage guarantee

Under exchangeability,

$$
\mathbb{P}\big(Y_{n+1} \in C(X_{n+1})\big) \ \ge\ 1 - \alpha
$$

exactly, for finite $n$. The $\lceil (n_c+1)(1-\alpha)\rceil$ correction is what
delivers finiteness: using the plain empirical quantile gives coverage only
asymptotically and **under-covers** at the $n \approx 20$–$50$ that
immuno-oncology cohorts actually have. With $n_c$ in the tens this is a
several-percentage-point difference, i.e. larger than most reported effect sizes.

This is a **marginal** guarantee. It does not promise coverage for any particular
subgroup — which matters, because the subgroup most likely to be under-covered is
the clinically decisive one. That is what §6.5 addresses.

### 6.5 Mondrian (group-conditional) conformal

Calibrate a separate quantile per group $g$ (here: per archetype):

$$
\hat q_g = s^{g}_{(\lceil (n_g+1)(1-\alpha)\rceil)}, \qquad
\text{score measured on calibration points whose true label is } g
$$

giving the strictly stronger guarantee

$$
\mathbb{P}\big(Y \in C(X) \ \big|\ Y = g\big) \ \ge\ 1-\alpha \quad \forall g
$$

Implementation note: at prediction time the true group is unknown, so each
candidate label $y$ is admitted against **its own** group quantile $\hat q_y$.
This makes an empty prediction set possible; an empty set is a coverage
violation, so the arg-max label is restored and the event is counted and
returned as an attribute rather than hidden.

Mondrian is on by default here. Marginal coverage can be satisfied while
systematically failing on the archetype that decides the clinical question, and a
tool that permits that failure by default would be reproducing the problem it
exists to detect.

### 6.6 Weighted conformal for covariate shift

When calibration and test distributions differ, exchangeability fails. With
importance weights $w(x) = p_{\text{test}}(x)/p_{\text{cal}}(x)$, replace the
empirical quantile by the **weighted** quantile:

$$
\hat q = \min\Big\lbrace q : \frac{\sum_{i:\ s_i \le q} w_i}{\sum_j w_j} \ \ge\ \frac{\lceil (n_c+1)(1-\alpha)\rceil}{n_c+1} \Big\rbrace 
$$

This restores approximate validity under the shift, and — importantly — the
*direction* of failure is informative: when the weights concentrate on a few
calibration points, $\hat q$ inflates and the model abstains more. That is the
correct behaviour, and it is why the conformal literature reports that coverage
degrades gracefully rather than silently.

### 6.7 Abstention rule and the selective-risk view

$$
\text{abstain}(x) = \mathbb{1}\big[\ |C(x)| > 1\ \big]
$$

Coverage alone is not a sufficient readout: a model can hold marginal coverage
while being non-committal for most samples. The package therefore reports the
**coverage–accuracy–abstention triple**, and draws the whole operating curve
(`plot_coverage_risk()`) rather than a single point:

- $\text{answered}(\tau)$ — fraction not abstained at threshold $\tau$;
- $\text{accuracy}(\tau)$ — accuracy **among those answered**;
- $\text{singleton rate}(\tau)$ — fraction on which the model committed.

### 6.8 Where recoverability deliberately is *not* invented

No synthetic "recoverability probability" is produced. A committed sample carries
its classifier confidence; an abstained sample carries `NA`, because no
calibrated call was made. Fabricating a score for an abstained sample would
reintroduce exactly the over-confidence the method exists to remove. The
empirical accuracy of the committed subset is reported separately
(`reliability_curve()`), so the number the user reads is measured, not asserted.

### 6.9 Reliability and expected calibration error

Bin calibration points by confidence. For bin $b$ with $n_b$ points:

$$
\text{ECE} = \sum_b \frac{n_b}{n}\ \big|\ \text{acc}_b - \overline{\text{conf}}_b\ \big|
$$

Plotted against the identity line: deviation above indicates over-confidence.

### 6.10 Pairwise archetype separability

For archetypes $i \ne j$:

$$
\Sigma_{ij} = \frac{1}{n_c}\sum_{t=1}^{n_c} \mathbb{1}\big[\lbrace i,j\rbrace \subseteq C_t\big]
$$

the fraction of calibration samples whose prediction set contained **both**
labels. High $\Sigma_{ij}$ means bulk data cannot separate that pair. This is the
actionable output: it names the clinical distinction the assay will not deliver,
instead of reporting one aggregate accuracy that hides it.

---

## 7. Sampling adequacy

### 7.1 Pseudo-bulk construction

For a region $R$ of diameter $d$ (a disc centred at a valid location containing
$m_R$ spots):

$$
x_R = \sum_{s \in R} c_s \qquad (\texttt{aggregate = "sum"})
\qquad\text{or}\qquad
x_R = \frac{1}{m_R}\sum_{s \in R} c_s
$$

with $c_s \in \mathbb{R}^p$ the spot's count vector. Sweeping $d$ generates
paired (bulk, spatial) data without new specimens — which is what makes
recoverability trainable at all. The region label is the **majority** spot label.

### 7.2 Hill saturation model

Performance against diameter follows a saturating curve, fitted by nonlinear
least squares:

$$
P(d) = P_\infty \frac{d^{\ h}}{K^{\ h} + d^{\ h}}
$$

$P_\infty$ is the plateau (asymptotic recoverable signal), $K$ the
half-saturation diameter, $h$ the cooperativity. A Hill form is used rather than a
single exponential because it accommodates a genuine **threshold** (large $h$) as
well as a gradual rise ($h \approx 1$); a one-parameter exponential forces a
shape and therefore biases the adequacy estimate.

### 7.3 Adequacy threshold

Solving $P(d^\star) = t\ P_\infty$ for target fraction $t$:

$$
d^\star = K \left(\frac{t}{1-t}\right)^{1/h}
$$

$d^\star$ is the **minimum region diameter** at which the spatial endpoint is
recoverable at fraction $t$ of its asymptote. At $t = 0.95$ this reproduces, as a
continuous and data-driven rule, the empirical observation from a randomised
neoadjuvant trial that sub-3 mm regions fail to preserve the predictive spatial
signal.

### 7.4 Why the scorer is deliberately simple

The default scorer is leave-one-out nearest-centroid accuracy. A heavier learner
would confound the sampling question with model capacity: the object of study is
how much signal *survives the sampling geometry*, not how well a flexible model
can eke out what remains. Model capacity is varied separately (§2–3), not here.

---

## 8. Computational summary

| Procedure | Cost | Dominant term |
|---|---|---|
| $\mathcal{V}(S)$ for one $S$ | $m$ fits + inner selection | $m \times G \times C_{\text{inner}}$ |
| Shapley, exact | $2^K$ cached evaluations | $2^K \times C_{\text{cv}}$ |
| Permutation null, one view | $2B$ evaluations | $2B \times C_{\text{cv}}$ |
| Full `attrib_views` | $1 + \tfrac{K B}{2^{K-1}}$ lattice passes | dominated by the null at small $K$ |
| Conformal fit | 1 fit + $O(n_c \log n_c)$ | the sort |
| Prediction | 1 fit + $O(n \log K)$ | — |
| Sampling sweep | $|D| \times |G|$ scorer evaluations | scorer is $O(n^2 p)$ |

$m$ = folds, $G$ = size of the $\alpha$ grid, $B$ = permutations, $D$ = diameters,
$G$ = panel sizes.

**Optimisations that are in place and were not obvious:**

1. **One lattice pass, reused across players** — the naive implementation
   recomputes $\mathcal{V}$ per player, costing $K$ times more.
2. **Shared folds across subsets** — removes both split noise and redundant fold construction.
3. **Canonical subset keys in view order** — the failure mode of getting this wrong is `subscript out of bounds`, not a wrong number, which is why the lattice is size-asserted.
4. **Frozen standardisation statistics** computed once per fold and applied to train and test — a correctness requirement, not merely a speed one.
5. **Deviance requested explicitly from `cv.glmnet`** — removes an $n$-dependent code path that silently changes what `cvm` means.

---

## 9. What these algorithms do **not** claim

Stated plainly, because over-claiming is the failure mode this package exists to
oppose:

- **Shapley attribution is associative, not causal.** A high $\phi_k$ means the
  layer is *predictively useful*, never that it drives the biology.
- **Conformal coverage is marginal unless Mondrian is enabled**, and degrades
  under distribution shift (§6.6). It is a calibrated *inference* statement, not
  a clinical guarantee.
- **Abstention is not accuracy.** A model that abstains on hard cases is honest,
  not better. `plot_coverage_risk()` exists so this cannot be read as a gain.
- **$d^\star$ is estimated from the data at hand.** It is a design guide derived
  from a specific tissue type and platform, and it should be re-estimated for a
  new platform rather than reused as a constant.
- **Nothing here substitutes for external validation.** Every quantity is
  estimated within the cohorts supplied.

---

## 10. References for the underlying methods

- Shapley LS (1953). *A value for n-person games.* Contributions to the Theory of Games II.
- Lindeman RH, Merenda PF, Gold RZ (1980). *Introduction to Bivariate and Multivariate Analysis* — the LMG regression decomposition.
- GrÖmping U (2006). Relative importance for linear regression in R: the package relaimpo. *J Stat Softw* 17(1).
- Vovk V, Gammerman A, Shafer G (2005). *Algorithmic Learning in a Random World* — conformal prediction.
- Lei J, G'Sell M, Rinaldo A, Tibshirani RJ, Wasserman L (2018). Distribution-free predictive inference for regression. *JASA* 113(523).
- Romano Y, Sesia M, Candès EJ (2020). Classification with valid and adaptive coverage. *NeurIPS* — APS.
- Boström H, Linusson H, Löfström T, Johansson U (2017). Accelerating difficulty estimation for conformal regression forests — Mondrian/group-conditional calibration.
- Tibshirani RJ, Barber RF, Candès EJ, Ramdas A (2019). Conformal prediction under covariate shift. *NeurIPS*.
- Hoffman GE, Schadt EE (2016). variancePartition. *Nucleic Acids Res* 44(18) — the tool this package's §3 is deliberately distinct from.
- Argelaguet R, et al. (2018, 2020). MOFA / MOFA+. *Mol Syst Biol* / *Genome Biol* — likewise.
