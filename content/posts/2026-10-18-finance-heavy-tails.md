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

### Decision: maximize mean return under an ES budget

At confidence $\alpha$, expected shortfall averages the worst $1-\alpha$ fraction of losses. For continuous losses,

$$
\operatorname{ES}_\alpha(L)=E[L\mid L\ge q_\alpha(L)],\qquad
\boxed{\max_{w\ge0,\,\mathbf1^\top w=1}\widehat\mu^\top w
\quad\text{s.t.}\quad\widehat{\operatorname{ES}}_\alpha(-w^\top r)\le b.}
$$

Here $\widehat\mu$ is the vector of training means; $b$ is the allowed **average tail loss**, rather than a maximum daily loss. At 97.5%, ES averages the worst 2.5% of outcomes.

The [Rockafellar–Uryasev (2000)](https://sites.math.washington.edu/~rtr/papers/rtr179-CVaR1.pdf) scenario formulation yields a linear program:

$$
\widehat{\operatorname{ES}}_\alpha(w)
=\min_\tau\left\{\tau+\frac{1}{n(1-\alpha)}
\sum_t(-w^\top r_t-\tau)_+\right\}.
$$

<details markdown="1">
<summary>Explicit optimization variables and constraints</summary>

Optimize $(w,\tau,s_1,\ldots,s_n)$ to maximize $\widehat\mu^\top w$, subject to

$$
\tau+\frac{\sum_ts_t}{n(1-\alpha)}\le b,\quad
s_t\ge-w^\top r_t-\tau,\quad s_t\ge0,\quad
w\ge0,\quad\mathbf1^\top w=1.
$$

The auxiliary $s_t$ records loss exceeding threshold $\tau$. The empirical calculation includes a fractional boundary observation when $n(1-\alpha)$ is noninteger.

</details>

### What changes if returns are assumed normal?

$$
\operatorname{ES}^{N}_\alpha(w)
=-w^\top\widehat\mu+
\sqrt{w^\top\widehat\Sigma w}\,
\frac{\phi\{\Phi^{-1}(\alpha)\}}{1-\alpha}.
$$

Historical ES uses **observed joint loss scenarios**. Normal ES uses **means and covariances**. With the same mean objective and budget, a smaller modelled tail risk permits more exposure to high-mean assets. The two allocations can therefore differ even on the same training data.

### Interactive lab B · allocation and frontier

**Start with “Equal weight,” then select “Optimize with historical ES.”** Move the budget to see which allocations become feasible. Switch to normal ES to compare the risk estimate with the observed tail loss of the same weights.

<iframe id="portfolio-lab" data-finance-lab src="{{base}}/assets/finance/portfolio.html" title="Ten-industry mean–ES allocation and efficient frontier" style="width:100%;height:1510px;border:1px solid #d8e2e9;border-radius:6px;" loading="lazy"></iframe>
<script src="{{base}}/assets/finance/embed.js" defer></script>

[Open portfolio lab full-width]({{base}}/assets/finance/portfolio.html)

### Read Figure 3 · capital weight differs from risk contribution

For the selected portfolio's tail event $\mathcal T_w=\{L_p\ge q_\alpha(L_p)\}$,

$$
\underbrace{C_j(w)}_{\text{industry }j\text{'s tail contribution}}
=w_jE[-r_j\mid\mathcal T_w],\qquad
\sum_{j=1}^{10}C_j(w)=\operatorname{ES}_\alpha(L_p).
$$

The left bars show capital weights; the right bars show each industry's losses **on the portfolio's bad days**. Changing weights also changes which days enter the tail. This is how joint downside behaviour affects diversification: an asset that gains on those days can have a negative contribution.

### Read Figure 4 · an efficient frontier

The feasible set contains all attainable pairs $(\operatorname{ES}_\alpha(w),\widehat\mu^\top w)$. Its efficient boundary retains portfolios for which no other allocation has both **lower or equal ES and higher or equal mean**, with one strict improvement.

$$
f(b)=\max_w\{\widehat\mu^\top w:\operatorname{ES}_\alpha(w)\le b,\ w\ge0,\ \mathbf1^\top w=1\}.
$$

Each displayed solution maximizes mean under one budget. The blue and dashed frontiers use different risk models; **the x-axis is ES, rather than standard deviation**. For normal-optimized weights, the horizontal orange segment connects their normal risk estimate to their historical ES. Both endpoints have the same mean and weights.

**What the exercise establishes:** tail modelling changes risk estimates, feasible allocations and the estimated frontier. The evaluation table shows what those fixed weights subsequently experienced. Heavy tails supply risk information; higher future returns depend on expected-return estimation and realised outcomes.

<details markdown="1">
<summary>Data, scope and reproducibility</summary>

- **Source:** [French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html), 6,539 daily observations in 2000–2025. Downloaded October 18, 2026; August 2026 CRSP vintage. Copyright Eugene F. Fama and Kenneth R. French.
- **Units:** equations are scale-consistent; all panel returns, means and ES values are daily percentages. Contributions are percentage points of loss. Industry series are research portfolios.
- **Estimation:** frontiers use 2010–2019 only; evaluation uses 2020–2025 with frozen target weights. Costs are excluded. Estimated means and tail expectations are uncertain, and future ES can exceed the training budget.
- **Tail interpretation:** these are unconditional empirical diagnostics. Volatility clustering can generate fat pooled tails. Simple-return losses on an unlevered long position are bounded at 100%; a polynomial tail is an idealized model over the relevant range.
- **Reproduce:** download [portfolio_study.py]({{base}}/assets/finance/portfolio_study.py) and [returns.json]({{base}}/assets/finance/returns.json) into one folder; install NumPy and SciPy, then run `python portfolio_study.py`. [Calculated weights and results]({{base}}/assets/finance/portfolio-results.json).

</details>
