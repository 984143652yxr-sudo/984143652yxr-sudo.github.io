---
{
  "slug": "finance-heavy-tails",
  "title": "Heavy tails in finance: how tail risk changes portfolio weights",
  "date": "2026-10-18",
  "date_label": "October 18, 2026",
  "status": "Interactive finance study",
  "reading_time": "10 min + exploration",
  "summary": "From one asset's extreme returns to portfolio tail contributions and maximum expected return under an expected-shortfall budget."
}
---

## 1. Whose distribution has a heavy tail?

**Heavy tails are a property of a return distribution.** They can appear at each level:

| Random variable | Interpretation |
|---|---|
| $r_j$ | Return of one stock |
| $r_A=\sum_jv_jr_j$ | Return of an industry portfolio |
| $r_p=w r_A+(1-w)r_B$ | Return of a portfolio of two industries |

The data here contain **industry portfolios**. Each is already a weighted collection of stocks. Two industries make allocation a one-dimensional problem: choose $w\in[0,1]$. The same model extends to many individual stocks.

For standardized returns $Z=(r-E[r])/\operatorname{SD}(r)$, compare

$$
\underbrace{P(|Z|>x)}_{\text{observed tail frequency}}
\quad\text{with}\quad
\underbrace{2[1-\Phi(x)]}_{\text{normal reference}}.
$$

Heavy tails relative to a normal benchmark mean unusually frequent large deviations. They do not determine the sign or size of $E[r]$. Portfolio allocation needs **expected returns, marginal tails, and joint dependence**.

A polynomial model has $P(X>x)\sim\ell(x)x^{-\kappa}$; smaller $\kappa$ means a heavier tail. Student's $t_\nu$ is a familiar example with exponent $\nu$. Empirical fat tails do not establish a power law or infinite variance. Long-only simple-return losses are bounded by 100%, so an unbounded tail model is an approximation over a relevant range. [Cont (2001)](https://rama.cont.perso.math.cnrs.fr/pdf/empirical.pdf)

## 2. Interactive lab: diagnose → combine → allocate

**First set A to 100%:** Figures 1–2 describe one industry's returns. Then mix A and B: Figures 3–4 show how their joint losses affect allocation. **Load worked example** sets technology–utilities, 2010–2019, 97.5% ES, and a 2.5% loss budget.

<iframe id="finance-lab" src="{{base}}/assets/finance/panel.html" title="Interactive heavy-tail and portfolio allocation laboratory" style="width:100%;height:1650px;border:1px solid #d8e2e9;border-radius:6px;" loading="lazy"></iframe>
<script src="{{base}}/assets/finance/embed.js" defer></script>

[Open full-width]({{base}}/assets/finance/panel.html) · [Download data]({{base}}/assets/finance/returns.csv)

| Figure | Axes | Statistical question |
|---|---|---|
| **1. Q–Q** | Normal quantile → observed standardized return | Are extreme returns larger than a normal model predicts? |
| **2. Exceedance** | Magnitude in SD units → probability, log scale | How much more frequent are large moves? |
| **3. Tail contributions** | A's loss → B's loss; one point per date | What does each asset do when **this portfolio** loses most? |
| **4. Allocation** | ES loss → estimated mean return; one point per weight | Which weight gives the highest estimated return within the loss budget? |

For the full-sample 50/50 technology–utilities portfolio, 100 of 6,539 days exceed three sample SDs in absolute value; the fitted normal reference gives about 17.7. This is an unconditional diagnostic, pooling calm and volatile periods.

## 3. Figure 3: marginal tails become portfolio risk through dependence

Let $L_A=-r_A$, $L_B=-r_B$, and $L_p=wL_A+(1-w)L_B$. Positive values are losses; negative values are gains.

Define the portfolio-tail event

$$
\mathcal T_w=\{L_p\ge\operatorname{VaR}_\alpha(L_p)\}.
$$

For a continuous loss distribution,

$$
\boxed{
\operatorname{ES}_\alpha(L_p)
=\underbrace{wE[L_A\mid\mathcal T_w]}_{C_A}
+\underbrace{(1-w)E[L_B\mid\mathcal T_w]}_{C_B}.
}
$$

**Orange points are the portfolio's worst days.** The panel reports $C_A$ and $C_B$ in percentage points of loss. It uses fractional weighting at the empirical tail boundary so their sum equals sample ES exactly.

- Large losses in both assets on orange days: limited protection from mixing them.
- B gains on those days: B can offset A's losses; $C_B$ can be negative.
- Change $w$: both the loss contributions and the dates in $\mathcal T_w$ can change.

Under smoothness, the allocation sensitivity is

$$
\frac{d\operatorname{ES}_\alpha(L_p)}{dw}
=E[L_A-L_B\mid\mathcal T_w].
$$

Thus the useful diversification question is **which asset reduces losses in the portfolio tail**, alongside its expected return. Two separate tail histograms cannot answer this joint question. [Dependence and risk](https://people.math.ethz.ch/~embrecht/ftp/pitfalls.pdf)

## 4. Figure 4: maximize return subject to tail risk

### Specify the decision

Let $\mu_A=E[r_A]$, $\mu_B=E[r_B]$, and $b$ be the allowed expected shortfall. The allocation problem is

$$
\boxed{
\max_{0\le w\le1}\;w\mu_A+(1-w)\mu_B
\quad\text{subject to}\quad
\operatorname{ES}_\alpha(wL_A+(1-w)L_B)\le b.
}
$$

Without a risk constraint, this linear objective puts everything in whichever asset has the larger mean. Tail risk restricts that concentration. A different objective, such as expected log wealth, can give different weights.

In Figure 4, choose the **highest point to the left of the budget line**. The star marks the historical-ES solution; the open diamond marks the normal-model solution. The button applies the historical solution to all four figures. An infeasible budget is reported explicitly.

### What changes when tails are retained?

For normally distributed portfolio returns,

$$
\operatorname{ES}^{N}_\alpha(L_p)
=-\mu_p+\sigma_p\frac{\phi(\Phi^{-1}(\alpha))}{1-\alpha}.
$$

This reduces the risk calculation to mean and variance. Historical ES instead retains observed tail magnitudes and same-day dependence. Its general definition and scenario representation are

$$
\operatorname{ES}_\alpha(L)
=\frac{1}{1-\alpha}\int_\alpha^1F_L^{-1}(p)\,dp,
$$

$$
\widehat{\operatorname{ES}}_\alpha(L_p)
=\min_\tau\left[\tau+
\frac{1}{n(1-\alpha)}\sum_t(L_{p,t}-\tau)_+\right].
$$

[Rockafellar & Uryasev (2000)](https://sites.math.washington.edu/~rtr/papers/rtr179-CVaR1.pdf) provide the classical CVaR optimization framework. For many assets, maximize $\widehat\mu^\top w$ subject to

$$
\tau+\frac{\sum_ts_t}{n(1-\alpha)}\le b,\quad
s_t\ge-w^\top r_t-\tau,\quad s_t\ge0,\quad
w\ge0,\quad\mathbf1^\top w=1.
$$

This is a linear program. ES incorporates extreme-loss severity without requiring a fitted power law. The panel solves the two-asset version on a 1% weight grid. All return and risk values shown in the panel are daily percentages.

### Worked example: the tail assumption changes the feasible weight

**Estimate on 2010–2019; freeze weights; evaluate on 2020–2025.** Both methods maximize the same estimated mean return with $\alpha=97.5\%$ and $b=2.5\%$.

| Risk model used to select weights | Technology | Utilities | Model training ES | Observed training ES | Test ES | Test mean return |
|---|---:|---:|---:|---:|---:|---:|
| Normal | 97% | 3% | 2.497% | 3.222% | 4.809% | 0.094% |
| Historical | 51% | 49% | 2.496% | 2.496% | 4.047% | 0.070% |

The normal model permits substantially more technology exposure. Its chosen allocation already exceeds the budget when measured by historical training ES. Retaining observed tail losses leads to a more balanced allocation, with lower test ES **and lower test mean return**.

Both test ES values exceed 2.5%: an estimated training constraint is not a guarantee about future risk. Results depend on the estimation window; using 2000–2019 with a 3.5% budget makes both methods select utilities alone. This worked example illustrates a mechanism, not universal superiority. [Reproduce this comparison]({{base}}/assets/finance/allocation_study.py) · [Exact results]({{base}}/assets/finance/allocation-results.json)

## 5. A second use: size total exposure

For risky-portfolio loss $L$, risky fraction $a\ge0$, and zero-return cash,

$$
\operatorname{ES}_\alpha(aL)=a\operatorname{ES}_\alpha(L),\qquad
 a\le\min\left(1,\frac{b}{\operatorname{ES}_\alpha(L)}\right).
$$

If estimated ES is 7% and the budget is 2%, exposure is at most $2/7\approx28.6\%$. A heavier fitted tail increases estimated ES and reduces allowed exposure.

A classical extension estimates time-varying volatility and then fits extreme-value tails to standardized residuals:

$$
r_t=\mu_t+\sigma_t\varepsilon_t
\quad\Longrightarrow\quad
\operatorname{ES}_\alpha(L_{t+1}\mid\mathcal F_t)
=-\mu_{t+1}+\sigma_{t+1}\operatorname{ES}_\alpha(-\varepsilon).
$$

[McNeil & Frey (2000)](https://statmath.wu.ac.at/~frey/publications/evt-garch.pdf) develop this GARCH–EVT approach. It provides conditional risk estimates for exposure control; the current panel uses historical unconditional ES.

## Data and estimation limits

- **Data:** ten U.S. value-weighted industry portfolios, 2000–2025, from the [French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html). Snapshot downloaded October 18, 2026; August 2026 CRSP vintage. Copyright Eugene F. Fama and Kenneth R. French.
- **Portfolio convention:** daily fixed weights, no leverage or transaction costs. Industry portfolios are research series, not ETFs. Equations use return fractions; tables display percentages.
- **Estimation:** sample means are noisy; tail samples are small. Figure 4 is in sample. Only the worked-example table uses frozen weights on later data. Use rolling validation and dependence-aware uncertainty estimates before drawing performance conclusions.
- **Code:** [data extraction]({{base}}/assets/finance/prepare_data.py), [allocation comparison]({{base}}/assets/finance/allocation_study.py), [interactive calculations]({{base}}/assets/finance/panel.js). Download `returns.json` and `prepare_data.py` beside `allocation_study.py` to reproduce the comparison.
