---
{
  "slug": "finance-heavy-tails",
  "title": "Heavy tails in finance: from return distributions to portfolio decisions",
  "date": "2026-10-18",
  "date_label": "October 18, 2026",
  "status": "Interactive finance study",
  "reading_time": "15 min + exploration",
  "summary": "Explore historical stock-industry returns, see departures from a normal model, and work through expected-shortfall allocation and extreme-value risk sizing."
}
---

## 1. Start with returns, then look at the tails

**Study question:** How does a model of rare, large losses change portfolio construction? A price chart shows a market's path. A return distribution asks how large a one-period change can be, and how often it occurs.

For a dividend-adjusted investment value $P_t$, the daily simple return is

$$
r_t=\frac{P_t-P_{t-1}}{P_{t-1}},\qquad
r_{p,t}=w r_{A,t}+(1-w)r_{B,t},\qquad L_t=-r_{p,t}.
$$

A positive $L_t$ is a loss. The portfolio equation assumes weights $w$ and $1-w$ at the start of each day. Holding those weights fixed across dates means daily rebalancing. The study excludes costs and leverage.

| Notation | Meaning |
|---|---|
| $r_t$, $r_{p,t}$ | Asset return and portfolio return |
| $L=-r_p$ | Portfolio loss as a fraction of initial wealth |
| $\mu,s$ | Mean and standard deviation of returns |
| $z_t=(r_{p,t}-\bar r_p)/s$ | Standardized observed return |
| $\alpha$ | Confidence level for a loss quantile or expected shortfall |
| $u,\xi,\beta$ | Extreme-value threshold, tail shape, and tail scale |

Equations use return fractions: $0.02$ means 2%. The interactive panel displays returns and losses in **percent**.

### Interactive lab

Start with the default 50/50 technology–utilities portfolio. Change an industry or allocation, compare a quiet year with a crisis year, and hover over extreme observations. The four plots update together.

<iframe id="finance-lab" src="{{base}}/assets/finance/panel.html" title="Interactive heavy-tail and portfolio allocation laboratory" style="width:100%;height:1550px;border:1px solid #d8e2e9;border-radius:6px;" loading="lazy"></iframe>
<script src="{{base}}/assets/finance/embed.js" defer></script>

[Open the lab in a full-width tab]({{base}}/assets/finance/panel.html) · [Download the historical returns]({{base}}/assets/finance/returns.csv)

**Data.** These are ten value-weighted U.S. stock-industry research portfolios from the [Kenneth R. French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html), including technology, energy, healthcare, and utilities. The fixed sample contains **6,539 trading days, January 3, 2000–December 31, 2025**. Returns include distributions; these are industry portfolios rather than individual stocks or investable ETFs. The downloaded file uses the August 2026 CRSP database. Historical revisions can change the results of a later download.

### How to read the plots

1. **Normal Q–Q:** the x-axis is a standard-normal quantile; the y-axis is an observed standardized return quantile. Under approximate normality, points follow the diagonal. Large downward and upward departures reveal unusually long lower and upper tails relative to that benchmark.
2. **Tail probability:** x is a magnitude threshold in standard-deviation units; y is the fraction of returns beyond that threshold in either direction. The logarithmic y-axis makes rare events visible. The gray line is $2[1-\Phi(x)]$.
3. **Joint returns:** x and y are the two industries' same-day returns. Orange points lie below both empirical fifth percentiles. This shows whether diversification partners experience large losses together.
4. **Allocation curve:** x is the weight in industry A; y is daily expected shortfall. The historical curve retains observed joint losses; the normal curve uses only the portfolio's estimated mean and standard deviation.

For the default full-sample portfolio, **100 days** exceed three sample standard deviations in absolute magnitude. A normal reference gives about **17.7 days**. This comparison is a descriptive diagnostic: parameters were estimated from the same sample, and observations are temporally dependent.

## 2. What does “heavy-tailed” mean?

In empirical finance, *fat tails* commonly describes more extreme returns than a fitted Gaussian model predicts. The [stylized-facts literature](https://rama.cont.perso.math.cnrs.fr/pdf/empirical.pdf) also emphasizes volatility clustering: large moves tend to occur near other large moves. A large sample standard deviation alone is different from an unusually heavy tail after standardization.

A more specific mathematical model is a regularly varying upper tail:

$$
P(X>x)\sim \ell(x)x^{-\kappa},\qquad \kappa>0,
$$

where $\ell$ is slowly varying. A smaller exponent $\kappa$ means a heavier tail. For example, Student's $t_\nu$ has a polynomial tail with exponent $\nu$, while the normal tail decreases approximately as $e^{-x^2/2}/x$. For a variance-matched comparison, use $\sqrt{(\nu-2)/\nu}\,t_\nu$ when $\nu>2$. Mean, variance, and fourth moment require $\nu>1,2,4$, respectively.

Regular variation is one class of heavy tails. A finite return sample cannot prove this asymptotic form or establish infinite variance. Moreover, an unlevered long-only simple-return loss is bounded by 100%. Polynomial models for such losses describe an observed range; they cannot be literally valid to infinity. Log returns and leveraged positions have different bounds.

### Why time variation matters

A useful decomposition is

$$
r_t=\mu_t+\sigma_t\varepsilon_t.
$$

Pooling calm and volatile days creates a scale mixture. Even conditionally Gaussian innovations can produce an unconditional distribution with pronounced tails. Consequently, two questions should be separated: **How extreme is a move relative to long-run volatility? How extreme is it relative to today's volatility?** The panel answers the first using the selected sample. It does not fit a conditional volatility model.

Marginal tails and joint dependence are also distinct. A stock can have large individual moves without always falling when another stock falls. The scatterplot preserves observed same-day dependence, which is essential when adding portfolio returns. Correlation alone does not describe general joint-tail behavior. [Embrechts, McNeil & Straumann](https://people.math.ethz.ch/~embrecht/ftp/pitfalls.pdf)

## 3. Classical application I: allocate by expected shortfall

### The target

Value-at-risk is a loss quantile; expected shortfall averages the upper loss tail:

$$
\operatorname{VaR}_\alpha(L)=F_L^{-1}(\alpha),\qquad
\operatorname{ES}_\alpha(L)=\frac{1}{1-\alpha}\int_\alpha^1F_L^{-1}(p)\,dp.
$$

At 97.5%, ES describes the average loss in the worst 2.5% probability mass. For continuous distributions it equals $E[L\mid L\ge\operatorname{VaR}_\alpha]$. The quantile-integral definition also handles discrete samples and ties.

For example, two strategies can share a 97.5% VaR of 3%, yet have ES values of 4% and 9%. VaR locates the threshold; ES distinguishes the severity beyond it. For normal portfolio returns,

$$
\operatorname{ES}^{N}_\alpha=-\mu_p+
\sigma_p\frac{\phi(\Phi^{-1}(\alpha))}{1-\alpha}.
$$

The panel compares this expression with empirical ES. It averages exactly the worst $n(1-\alpha)$ observations' worth of probability, using a fractional boundary weight when necessary.

### From a risk measure to an optimization problem

[Rockafellar & Uryasev (2000)](https://sites.math.washington.edu/~rtr/papers/rtr179-CVaR1.pdf) gave a practical scenario-based formulation of CVaR/ES optimization. With historical return vectors $r_t$ and portfolio weights $w$, solve

$$
\min_{w,\tau}\left\{
\tau+\frac{1}{n(1-\alpha)}\sum_{t=1}^{n}(-w^\top r_t-\tau)_+
\right\},\qquad
w\ge0,\quad\mathbf1^\top w=1.
$$

Here $\tau$ is an auxiliary loss threshold, and $(x)_+=\max(x,0)$. For many assets, introduce $s_t\ge0$ and $s_t\ge-w^\top r_t-\tau$; the objective becomes $\tau+\sum_ts_t/[n(1-\alpha)]$, a linear program. A required mean return can be added as $w^\top\widehat\mu\ge m_0$.

This approach lets severe observed or simulated joint losses influence weights directly. It requires an integrable loss distribution but no particular power-law assumption. Large tail losses matter more than their frequency alone.

### A reproducible allocation example

Use technology and utilities, choose weights on **2000–2019 only**, and evaluate the fixed allocation on **2020–2025**. Compare equal weights, minimum sample variance, and minimum 97.5% ES on a 1-percentage-point weight grid. There is no target-return constraint here, so this is a risk-minimization exercise.

| Training rule | Technology weight | Training ES | Test daily SD | Test ES |
|---|---:|---:|---:|---:|
| Equal weights | 50% | 3.455% | 1.342% | 4.037% |
| Minimum variance | 19% | 3.280% | 1.306% | 3.988% |
| Minimum ES | 20% | 3.280% | 1.304% | 3.980% |

The two optimized allocations are nearly identical. Their test ES difference is about **0.008 percentage points**, with no significance claim. This example shows how to formulate and evaluate a tail-risk objective; it does not establish superior future performance. Both optimized allocations also have lower test average daily returns than the 50/50 allocation (about 0.053–0.054% versus 0.069%).

The interactive allocation curve re-estimates risk within whichever period is selected. Its minimum is **in sample**. The table above uses a fixed chronological split and contains 5,031 training days and 1,508 test days. [Computation and source code]({{base}}/assets/finance/prepare_data.py) · [Exact results]({{base}}/assets/finance/analysis.json)

## 4. Classical application II: extreme-value tails for position sizing

Historical ES cannot generate a loss more severe than the observed scenarios. Extreme-value theory provides a model for the part of a loss distribution beyond a high threshold.

### Model the exceedances

For losses $L>u$, approximate the excess $Y=L-u$ by a generalized Pareto distribution (GPD):

$$
P(Y>y\mid L>u)\approx
\left(1+\xi\frac{y}{\beta}\right)^{-1/\xi},\qquad y\ge0,
$$

on its valid support, with scale $\beta>0$. Positive $\xi$ gives a polynomial tail; $\xi=0$ is the exponential limit. Let $p_u=P(L>u)$ be estimated by the exceedance fraction. For $\alpha>1-p_u$ and $\xi\ne0$,

$$
\operatorname{VaR}_\alpha\approx
u+\frac{\beta}{\xi}\left[\left(\frac{p_u}{1-\alpha}\right)^\xi-1\right],
\qquad
\operatorname{ES}_\alpha\approx
\frac{\operatorname{VaR}_\alpha+\beta-\xi u}{1-\xi},\quad\xi<1.
$$

ES is finite under the fitted GPD only when $\xi<1$.

### Turn a tail estimate into an exposure limit

Consider an **illustrative model**, separate from the historical data above: $u=3\%$, $p_u=5\%$, and $\beta=1\%$. Vary the tail shape while keeping the threshold, exceedance probability, and scale fixed.

| Tail shape $\xi$ | 99% VaR | 99% ES | Maximum risky fraction for a 2% ES budget |
|---|---:|---:|---:|
| 0.10 | 4.75% | 6.05% | 33.1% |
| 0.25 | 4.98% | 6.98% | 28.7% |
| 0.50 | 5.47% | 9.94% | 20.1% |

Place fraction $a$ in a given risky portfolio and the remainder in zero-return cash for the one-day horizon. Positive homogeneity gives

$$
\operatorname{ES}_\alpha(aL)=a\operatorname{ES}_\alpha(L),\qquad
 a\le\min\left(1,\frac{\text{loss budget}}{\operatorname{ES}_\alpha(L)}\right).
$$

Heavier modeled tails reduce the allowed exposure for a fixed risk budget. This controls the model's expected tail loss, rather than imposing a maximum realized loss. The numbers are hypothetical calculations, not estimated allocations or investment recommendations.

### A classical conditional version

[McNeil & Frey (2000)](https://statmath.wu.ac.at/~frey/publications/evt-garch.pdf) combine GARCH volatility filtering with extreme-value estimation of innovation tails. The analysis sequence is

$$
\text{return history}\longrightarrow(\widehat\mu_t,\widehat\sigma_t)
\longrightarrow\widehat\varepsilon_t
\longrightarrow\text{GPD on }-\widehat\varepsilon_t
\longrightarrow\text{next-day loss risk}.
$$

Given $r_{t+1}=\mu_{t+1}+\sigma_{t+1}\varepsilon_{t+1}$ and $X=-\varepsilon$,

$$
\operatorname{ES}_\alpha(L_{t+1}\mid\mathcal F_t)
=-\mu_{t+1}+\sigma_{t+1}\operatorname{ES}_\alpha(X).
$$

This separates the current volatility level from the innovation-tail shape. A risk budget can then be updated using the conditional estimate. The panel reports unconditional historical diagnostics; it does not claim to implement this GARCH–EVT fit.

Threshold selection trades approximation bias against the number of exceedances. Check stability across thresholds, residual dependence, and out-of-sample exceedances. Position limits also depend on liquidity, costs, and model uncertainty. A favorable tail fit on yesterday's data does not make extrapolation certain.

## 5. A small study workflow

1. Compare **2000–2025** with **2017**, **2008**, and **2020**. Record how the Q–Q plot and tail frequency change after standardization.
2. Compare technology–utilities with technology–energy. Use the joint-return plot to inspect simultaneous losses before interpreting diversification.
3. Set ES confidence to 95%, 97.5%, and 99%. Record the effective number of tail observations; high confidence levels can leave very little data in a single year.
4. Select weights on a training interval, freeze them, and assess losses on later dates. Use rolling tests and a block bootstrap for uncertainty, preserving serial and same-day cross-asset dependence.

The notebook provides historical diagnostics, one computed allocation exercise, and one hypothetical extreme-value sizing example. It makes no claim that heavy tails alone create a trading profit.

## Sources and reproducibility

- Kenneth R. French. [Ten Industry Portfolios: construction and coverage](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_10_ind_port.html). [Original daily CSV archive](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip). Copyright Eugene F. Fama and Kenneth R. French.
- Cont (2001). [Empirical properties of asset returns: stylized facts and statistical issues](https://rama.cont.perso.math.cnrs.fr/pdf/empirical.pdf).
- Rockafellar & Uryasev (2000). [Optimization of Conditional Value-at-Risk](https://sites.math.washington.edu/~rtr/papers/rtr179-CVaR1.pdf).
- McNeil & Frey (2000). [Estimation of tail-related risk measures for heteroscedastic financial time series: an extreme value approach](https://statmath.wu.ac.at/~frey/publications/evt-garch.pdf).
- Embrechts, McNeil & Straumann (2002). [Correlation and dependence in risk management: properties and pitfalls](https://people.math.ethz.ch/~embrecht/ftp/pitfalls.pdf).

Snapshot downloaded October 18, 2026. Missing-data codes are excluded; observations retain their trading-date alignment. The [data file]({{base}}/assets/finance/returns.json) records the original archive hash and version header. The [Python script]({{base}}/assets/finance/prepare_data.py) reproduces the data extraction and chronological allocation comparison; the [panel calculations]({{base}}/assets/finance/panel.js) are available for inspection. Returns are rounded to the precision supplied by the source.
