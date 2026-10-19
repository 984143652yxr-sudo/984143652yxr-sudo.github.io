---
{
  "slug": "finance-heavy-tails",
  "title": "Heavy tails: return diagnostics and portfolio allocation",
  "date": "2026-10-18",
  "date_label": "October 18, 2026",
  "status": "Interactive finance study",
  "reading_time": "7 min + exploration",
  "summary": "Diagnose heavy tails, then compare rolling ES and volatility budgets through next-month returns, cumulative wealth and drawdown."
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
| Data in this exercise | French value-weighted daily industry returns; estimate on rolling 36-month windows, evaluate on 2020–2025 |
| Allocation | Ten industry weights plus cash, $w_j\ge0$, $\sum_jw_j=1$; rebalance monthly |
| Relation to the paper | A new rolling risk-budget exercise on its 10Ind universe; the published experiment used monthly data and different estimation procedures |

### Rolling decision: maximize return under a risk budget

At each month-end $t$, estimate from the **previous 36 calendar months** $\mathcal D_t$. Choose weights, hold through the next month, then repeat. The evaluation period contains **72 months, January 2020–December 2025**.

$$
\mathcal D_t\ \longrightarrow\ (\widehat\mu_t,\widehat\Sigma_t)
\ \longrightarrow\ w_t\ \longrightarrow\ R_{t+1}^{\mathrm{out}}.
$$

The investment universe includes ten industries and **zero-return cash**, allowing low risk budgets to remain feasible. All strategies are long-only and unlevered, with $w\ge0$, $\mathbf1^\top w=1$. Set $\alpha=97.5\%$ and define daily risk in percentage units:

$$
m_t(w)=\widehat\mu_t^\top w,\qquad
\sigma_t(w)=\sqrt{w^\top\widehat\Sigma_t w},\qquad
N_t(w)=-m_t(w)+k_\alpha\sigma_t(w),\quad
k_\alpha=\frac{\phi(\Phi^{-1}(\alpha))}{1-\alpha}.
$$

Historical ES uses the [Rockafellar–Uryasev scenario formulation](https://sites.math.washington.edu/~rtr/papers/rtr179-CVaR1.pdf):

$$
H_t(w)=\min_\tau\left[\tau+
\frac{1}{n_t(1-\alpha)}\sum_{s\in\mathcal D_t}(-w^\top r_s-\tau)_+\right].
$$

| Strategy | Monthly decision | Information |
|---|---|---|
| **H · historical ES** | $w_t^H(b)=\arg\max_w m_t(w)$ subject to $H_t(w)\le b$ | Previous 36 months |
| **N · normal ES** | $w_t^N(b)=\arg\max_w m_t(w)$ subject to $N_t(w)\le b$ | Previous 36 months |
| **M · Markowitz** | $w_t^M(v)=\arg\max_w m_t(w)$ subject to $\sigma_t(w)\le v$ | Previous 36 months |
| **O · hindsight** | Maximize next-month terminal wealth subject to that month's realised ES $\le b$ | Entire next month, known in advance |

H and N share ES budget $b$. M uses a separate daily volatility budget $v$. Its initial setting $v=b/k_\alpha$ is a zero-mean normal scale conversion; equal numerical scale does not ensure equal realised tail risk. [Markowitz (1952)](https://doi.org/10.1111/j.1540-6261.1952.tb01525.x)

### Figure 4 · estimated risk becomes realised performance

**4A:** average the 72 rolling-window forecasts. Gray → orange keeps the normal-chosen weights and replaces their risk estimate. **4B:** apply those exact monthly weights to the subsequent month's actual returns, then pool all evaluation days. Move the ES budget to select corresponding blue and orange points in both panels; the Markowitz volatility slider is independent. Formulas for every line appear directly above its plot.

<iframe id="portfolio-lab" data-finance-lab src="{{base}}/assets/finance/portfolio.html?v=rolling36-1" title="Rolling Figure 4: estimated risk, realised returns, wealth and drawdown" style="width:100%;height:4400px;border:1px solid #d8e2e9;border-radius:6px;" loading="lazy"></iframe>
<script src="{{base}}/assets/finance/embed.js" defer></script>

[Open Figure 4 full-width]({{base}}/assets/finance/portfolio.html?v=rolling36-1)

### Monthly holdings and the future-informed optimum

Let $G_{t,d,j}$ be asset $j$'s cumulative gross return from the start of holding month $t+1$ to day $d$, with $G_{t,0,j}=1$. Opening weights $w_t$ remain invested; daily weights drift naturally. The actual daily portfolio return, in percent, is

$$
R_{t,d}(w_t)=100\left(\frac{G_{t,d}^\top w_t}{G_{t,d-1}^\top w_t}-1\right).
$$

The green benchmark observes the whole next month and solves

$$
w_t^O(b)=\arg\max_{w\ge0,\,\mathbf1^\top w=1}G_{t,\mathrm{end}}^\top w
\quad\text{s.t.}\quad
\widehat{ES}_{\alpha,t+1}\{-R_{t,d}(w)\}\le b.
$$

Each month has fewer than 40 trading days, so at $\alpha=97.5\%$ its empirical ES equals its largest daily loss. The constraint therefore becomes linear:

$$
\left[G_{t,d}-(1-b/100)G_{t,d-1}\right]^\top w\ge0
\quad\text{for every day }d.
$$

This yields an **exact monthly terminal-wealth optimum under monthly ES control**. Green uses future information and a terminal-wealth objective; it is a perfect-information reference. It does not estimate an achievable population optimum. The other three strategies use past observations only.

### Read three outcomes separately

1. **Risk calibration:** in 4A, compare normal and historical risk estimates for identical N weights. In 4B, inspect the later realised ES. Averaging window estimates and measuring ES over the entire backtest are different operations.
2. **Allocation performance:** compare H, N and M using both realised ES and realised mean. Lower estimated risk can permit larger positions; the actual tail loss depends on subsequent returns. The realised curves can bend and need not be efficient frontiers.
3. **Investment path:** cumulative wealth compounds the actual daily returns. Maximum drawdown measures the largest peak-to-trough loss, while ES averages the worst 2.5% of daily losses. These measure different aspects of risk.

<details markdown="1">
<summary>Data, assumptions and reproducibility</summary>

- **Source:** [French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html), value-weighted daily ten-industry returns. Downloaded October 18, 2026; August 2026 CRSP vintage. Copyright Eugene F. Fama and Kenneth R. French.
- **Timing:** January 2020 uses January 2017–December 2019 for estimation. Each subsequent window moves forward one month. Month-boundary trades are assumed frictionless; positions are held until the next rebalance.
- **Risk model:** training scenarios use linear one-day portfolio returns $w^\top r_s$. The evaluation implements monthly holdings with within-month drift. Historical ES and normal ES limits are fitted constraints; future realised losses can exceed them.
- **Cash and costs:** cash earns zero; no leverage, transaction costs, taxes or execution slippage. The industry series are research portfolios. The curves summarize 17 risk-budget settings.
- **Hindsight:** the green rule sees next-month prices and optimizes every month under its own realised ES constraint. At this confidence level, its daily losses are bounded by $b$. Its result is an optimistic reference with explicit future information.
- **Diagnostics:** the first laboratory describes unconditional sample tails; volatility clustering can contribute to heavy-looking pooled tails. An unlevered long position's simple loss is bounded at 100%.
- **Reproduce:** download [rolling_study.py]({{base}}/assets/finance/rolling_study.py) and [returns.json]({{base}}/assets/finance/returns.json), install NumPy and SciPy, then run `python rolling_study.py`. [Complete monthly weights, forecasts and paths]({{base}}/assets/finance/rolling-results.json).

</details>
