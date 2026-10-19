---
{
  "slug": "finance-heavy-tails",
  "title": "Heavy tails: return diagnostics and portfolio allocation",
  "date": "2026-10-18",
  "date_label": "October 18, 2026",
  "status": "Interactive finance study",
  "reading_time": "7 min + exploration",
  "summary": "Two separate laboratories: diagnose extreme returns, then allocate across a published ten-industry universe using an expected-shortfall frontier."
}
---

## 1. Heavy tails: what are they, and why do they matter?

**Intuition.** Heavy tails mean unusually frequent large returns, compared with a normal distribution of the same mean and variance.

$$
r_t=\frac{P_t-P_{t-1}}{P_{t-1}},\qquad
Z=\frac{r-\mu}{\sigma},\qquad
P(|Z|>x)>2\{1-\Phi(x)\}\quad\text{in the observed tail range}.
$$

**Model.** Student's $t_\nu$ has polynomial tails and finite variance for $\nu>2$; smaller $\nu$ means heavier tails.

$$
P(|T_\nu|>x)\sim Cx^{-\nu},\qquad
\log P(|T_\nu|>x)\sim-\nu\log x+\log C.
$$

**Diagnostic.** Compare empirical and normal exceedance probabilities; this describes the sample tail. [Cont (2001)](https://rama.cont.perso.math.cnrs.fr/pdf/empirical.pdf)

$$
\widehat P(|Z|>x)=\frac1n\sum_{t=1}^{n}\mathbf1\{|z_t|>x\},
\qquad P(|N(0,1)|>3)\approx0.27\%.
$$

### Interactive lab A · one distribution

Select an industry, then compare a quiet year with a crisis year. **Figure 1 measures the size of extremes; Figure 2 measures their frequency.** All controls here concern the return distribution.

<span id="interactive-lab"></span>
<iframe id="tail-lab" data-finance-lab src="{{base}}/assets/finance/panel.html" title="Heavy-tail diagnostics for one industry" style="width:100%;height:830px;border:1px solid #d8e2e9;border-radius:6px;" loading="lazy"></iframe>

[Open diagnostics full-width]({{base}}/assets/finance/panel.html)

### Why a strategy can fail during an extreme move

For wealth $W_t$ and risky exposure $a_t$ chosen before the next return,

$$
W_{t+1}=W_t(1+a_t r_{t+1}),\qquad
 a_t=\frac{\sigma_{\mathrm{target}}}{\widehat\sigma_t}
\quad\text{for a volatility-targeting rule with zero-return cash}.
$$

A small volatility estimate permits a large exposure. A tail shock arrives **before the next rebalance**: at $a_t=2$, a $-10\%$ asset return implies a $-20\%$ wealth return, before funding and trading costs. Matching volatility alone leaves the extreme-loss probability unspecified.

**Historical case: portfolio insurance in October 1987.** The S&P 500 fell about 20% on October 19. Dynamic insurance called for reducing equity exposure as prices fell; delayed execution and market illiquidity made the adjustment difficult. [Federal Reserve study](https://www.federalreserve.gov/pubs/feds/2007/200713/index.html)

<div style="display:flex;flex-wrap:wrap;gap:8px;margin:18px 0;font-size:14px;">
<div style="flex:1;min-width:135px;padding:12px;background:#f0f6f8;"><strong>1. Large price fall</strong><br>Existing exposure loses value</div>
<div style="flex:1;min-width:135px;padding:12px;background:#f0f6f8;"><strong>2. Protection rule</strong><br>Sell stocks or futures</div>
<div style="flex:1;min-width:135px;padding:12px;background:#fbf3eb;"><strong>3. Thin market</strong><br>Many sellers, delayed fills</div>
<div style="flex:1;min-width:135px;padding:12px;background:#fbf3eb;"><strong>4. Further pressure</strong><br>Protection depends on execution</div>
</div>

This illustrates two modelling needs: **tail scenarios for losses**, and **liquidity assumptions for rebalancing**. ES addresses loss severity; execution risk requires an additional model.

## 2. Portfolio allocation: use the tail information

### A published portfolio universe

Use the **ten U.S. industry portfolios (10Ind)** in Table 1 of [DeMiguel, Garlappi, Nogales & Uppal (2009)](https://pubsonline.informs.org/doi/10.1287/mnsc.1080.0986), with the familiar $1/N$ benchmark:

$$
w^{\mathrm{EW}}=(0.1,\ldots,0.1)^\top,\qquad
r_{p,t}=w^\top r_t,\qquad L_{p,t}=-w^\top r_t.
$$

Each component of $r_t\in\mathbb R^{10}$ is an industry's daily return. Each industry already holds many stocks. Positive $L$ denotes a loss.

| Component | Specification |
|---|---|
| Published universe | Consumer non-durables, durables, manufacturing, energy, technology, telecom, shops, healthcare, utilities, other |
| Data in this exercise | French value-weighted daily industry returns; estimate on 2010–2019, evaluate on 2020–2025 |
| Allocation | Ten weights, $w_j\ge0$, $\sum_jw_j=1$; maintain weights by daily rebalancing |
| Relation to the paper | A new mean–ES exercise on its 10Ind universe; the published experiment used monthly data and different estimation procedures |

### One target, two allocation rules

Let $m(w)=\widehat\mu^\top w$, $\sigma(w)=\sqrt{w^\top\widehat\Sigma w}$ and $\mathcal W=\{w\ge0:\mathbf1^\top w=1\}$. All estimates below use **2010–2019**. Confidence $\alpha=97.5\%$ selects the portfolio's worst 2.5% of days.

$$
\underbrace{H_\alpha(w)}_{\text{historical CVaR / ES}}
=\min_\tau\left[\tau+\frac{1}{n(1-\alpha)}\sum_t(-w^\top r_t-\tau)_+\right],\qquad
\underbrace{N_\alpha(w)}_{\text{normal ES}}=-m(w)+k_\alpha\sigma(w),\quad
k_\alpha=\frac{\phi(\Phi^{-1}(\alpha))}{1-\alpha}.
$$

The [Rockafellar–Uryasev formulation](https://sites.math.washington.edu/~rtr/papers/rtr179-CVaR1.pdf) minimizes the scenario tail average through a linear program. It computes empirical ES exactly, including a fractional boundary day; future risk still needs validation.

$$
\boxed{w_E(q)=\arg\min_{w\in\mathcal W,\,m(w)\ge q}H_\alpha(w)},\qquad
\boxed{w_M(q)=\arg\min_{w\in\mathcal W,\,m(w)\ge q}\sigma^2(w)}.
$$

**E** controls tail loss; **M** is [Markowitz minimum variance](https://doi.org/10.1111/j.1540-6261.1952.tb01525.x). The common return target $q$ makes their objectives comparable. At exactly fixed mean $m(w)=q$, minimizing normal ES also minimizes variance: the normal and Markowitz efficient families coincide. Markowitz optimization itself only requires finite covariance.

### Figure 4 · risk estimates, allocation quality and hindsight

**Figure 3** displays the selected weights. **Figure 4A** compares training risk; **Figure 4B** evaluates the same weights on later data. Both belong to the same interactive Figure 4 and share the controls.

<iframe id="portfolio-lab" data-finance-lab src="{{base}}/assets/finance/portfolio.html?v=combined-4" title="Figure 4: normal, CVaR, Markowitz and hindsight" style="width:100%;height:1900px;border:1px solid #d8e2e9;border-radius:6px;" loading="lazy"></iframe>
<script src="{{base}}/assets/finance/embed.js" defer></script>

[Open Figure 4 full-width]({{base}}/assets/finance/portfolio.html?v=combined-4)

Use $H'_\alpha,m'$ for **2020–2025** empirical ES and mean. The green hindsight curve solves

$$
w_O(u)=\arg\min_{w\in\mathcal W,\,m'(w)\ge u}H'_\alpha(w).
$$

It uses the entire later period to choose weights. This provides an optimistic empirical benchmark; the population optimum remains unknown.

| Line | Figure 4A coordinates | Figure 4B coordinates |
|---|---|---|
| Gray dashed · normal estimate | $(N_\alpha(w_M),m(w_M))$ | — |
| Orange dashed · Markowitz | $(H_\alpha(w_M),m(w_M))$ | $(H'_\alpha(w_M),m'(w_M))$ |
| Blue · CVaR | $(H_\alpha(w_E),m(w_E))$ | $(H'_\alpha(w_E),m'(w_E))$ |
| Green dotted · hindsight | — | $(H'_\alpha(w_O(u)),m'(w_O(u)))$ |

**1. Does normal overlook tail risk?** Follow gray A → orange M: identical weights, different estimates.

$$
\Delta_N=H_\alpha(w_M)-N_\alpha(w_M).
$$

A positive gap means normal **underestimates tail loss and overstates portfolio attractiveness**. Orange being close to blue concerns allocation efficiency; gray being far left concerns risk calibration. These are different effects.

**2. What does CVaR gain over variance control?** At the same training target,

$$
\Delta_{\mathrm{allocation}}=H_\alpha(w_M)-H_\alpha(w_E)\ge0.
$$

At the default target, training ES is **2.452% for CVaR versus 2.475% for Markowitz**. The gain is small. Later, Markowitz has lower ES (**3.605% versus 3.804%**) and lower mean (**0.0493% versus 0.0581% per day**). Tail optimization improves its training objective; future outcomes also depend on estimation and market changes.

**3. How close is the estimated tail risk to what occurred?** Keep E's weights fixed:

$$
\Delta_{\mathrm{calibration}}=H'_\alpha(w_E)-H_\alpha(w_E),\qquad
\Delta_{\mathrm{hindsight}}=H'_\alpha(w_E)-H'_\alpha\!\left(w_O(m'(w_E))\right).
$$

Default: training ES **2.452% → later ES 3.804%**, a **1.351 pp** increase. Hindsight achieves **3.655%** at at least E's later mean, a separate **0.148 pp** allocation improvement. The first compares estimated and realised risk; the second compares weights with and without future information. Later empirical ES is a finite-sample outcome, and regime changes contribute to its difference from training ES.

**Reading direction:** left is lower tail loss; up is higher mean. Compare allocations within one panel using historical ES. Each portfolio's tail is its own worst $1-\alpha$ fraction of days. Curves connect fitted targets; later paths may lose their training efficiency.

[Reproduction script]({{base}}/assets/finance/hindsight_study.py) · [Full weights and results]({{base}}/assets/finance/hindsight-results.json).

<details markdown="1">
<summary>Data, scope and reproducibility</summary>

- **Source:** [French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html), 6,539 daily observations in 2000–2025. Downloaded October 18, 2026; August 2026 CRSP vintage. Copyright Eugene F. Fama and Kenneth R. French.
- **Units:** equations are scale-consistent; all panel returns, means and ES values are daily percentages. Contributions are percentage points of loss. Industry series are research portfolios.
- **Estimation:** frontiers use 2010–2019 only; evaluation uses 2020–2025 with frozen target weights. Costs are excluded. Estimated means and tail expectations are uncertain, and future ES can exceed the training budget.
- **Tail interpretation:** these are unconditional empirical diagnostics. Volatility clustering can generate fat pooled tails. Simple-return losses on an unlevered long position are bounded at 100%; a polynomial tail is an idealized model over the relevant range.
- **Reproduce:** download [hindsight_study.py]({{base}}/assets/finance/hindsight_study.py) and [returns.json]({{base}}/assets/finance/returns.json) into one folder; install NumPy and SciPy, then run `python hindsight_study.py`. [Calculated weights and results]({{base}}/assets/finance/hindsight-results.json).

</details>
