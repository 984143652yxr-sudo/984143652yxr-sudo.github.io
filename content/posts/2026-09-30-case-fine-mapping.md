---
{
  "slug": "case-fine-mapping",
  "title": "From CIGMA to CASE: which SNP, in which cell type?",
  "date": "2026-09-30",
  "date_label": "September 30, 2026",
  "status": "Reading & local investigation",
  "reading_time": "13 min read",
  "summary": "Moving from gene-level genetic variance to SNP-by-cell-type fine-mapping: my CASE reading notes, an LD example, and a preliminary OneK1K overlap analysis."
}
---

## 1. The question left after CIGMA

After studying [CIGMA on OneK1K]({{base}}/diary/cigma-onek1k/), my next question is: **which variant contributes to a gene's cell-type-dependent regulation, and in which cells does it act?**

CIGMA quantifies shared and cell-type-specific genetic variance for a gene and can detect cell-type-specific eGenes. Its variance components summarize a set of SNPs. Individual-variant localization and GWAS colocalization require additional inference. This is the gap that motivates my study of **CASE: Cell-type-specific And Shared EQTL fine-mapping**. [CIGMA paper](https://www.nature.com/articles/s41586-026-10577-6)

For one gene, the simplified CIGMA construction from my [previous note]({{base}}/diary/cell-study/) is

$$
Y_{ic}=\mu_c+W_i^\top\gamma_c+\sum_j Z_{ij}(\alpha_j+\eta_{jc})+r_{ic}.
$$

The shared component $\alpha_j$ and type-dependent deviation $\eta_{jc}$ induce genetic variances $\sigma^2_{\rm shared}$ and $v_c$. Their aggregate specificity is

$$
\text{specificity}=\frac{\bar v}{\sigma^2_{\rm shared}+\bar v},
\qquad \bar v=\text{mean of }v_c\text{ across cell types}.
$$

This simplified equation keeps the estimand visible; the fitted CIGMA model additionally handles its specified residual and sampling structure.

A variance summary leaves many possible SNP-effect configurations compatible with the result. For example, one strongly varying SNP and several weakly varying SNPs can contribute similar total variance. The follow-up therefore needs a model for the **individual effect matrix**.

| Question | Output I need |
|---|---|
| How much of a gene's genetic regulation varies across cell types? | Shared/specific variance components |
| Which SNPs plausibly generate the signal? | Variant posterior probabilities and credible sets |
| In which cell types is each SNP supported? | A SNP-by-cell-type effect profile |
| Does the regulatory signal share a causal variant with a disease association? | A separate colocalization analysis using compatible locus-level data |

### Why this distinction matters quantitatively

In the published CIGMA OneK1K analysis, 193 of 10,288 genes passed the cs-eGene threshold, using 928 donors and seven cell types. Estimated specificity was 31.1% for cis and 59.4% for trans regulation. A bulk-like aggregate gave about 5% cis specificity, versus about 30% at cell-type resolution. These are variance fractions; the number of significant genes also depends on detection power. [CIGMA results](https://www.nature.com/articles/s41586-026-10577-6)

CASE reported 5,057 eGenes among 11,704 candidate genes in OneK1K: 13.0% more than mvSuSiE and 32.9% more than SuSiE. Those comparisons concern fine-mapping discoveries under the paper's design. [Lin et al., CASE (2026)](https://www.nature.com/articles/s41467-026-72176-3)

The two counts answer different questions: CIGMA's cs-eGene test concerns heterogeneity, while CASE's eGene call requires a credible set in at least one cell type. Their denominators and preprocessing also differ. My practical motivation is to connect **heterogeneous genetic regulation → candidate variants → cellular context**.

## 2. My notes: LD can imitate or hide sharing

The two sketches in my handwritten notes are the starting point. A marginal SNP association contains contributions from other correlated SNPs. Consequently, a pattern of marginal associations across cell types can differ from the pattern of underlying joint effects.

For standardized genotype columns and one cell type,

$$
y=X\beta+\epsilon,\qquad
\widehat b=X^\top y/N,
\qquad E[\widehat b\mid\beta]\approx R\beta,
$$

where $R=X^\top X/N$ is the LD correlation matrix under this normalization. Across cell types, this becomes $E[\widehat B\mid B]\approx RB$.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/ld-effects.png"><img src="{{base}}/assets/case-study/ld-effects.png" alt="Two deterministic examples comparing joint SNP effects with LD-mixed marginal associations across two cell types"></a>
<figcaption>My numerical version of the two LD sketches. Rows are SNPs, columns are cell types; entries are effects in arbitrary teaching units. The right panels are computed as RB. Click to enlarge.</figcaption>
</figure>

**Top row — apparent sharing.** SNP 1 acts only in cell type 1; SNP 2 acts only in cell type 2. With LD correlation 0.8, each SNP has a marginal association in both types. Interpreting those marginal associations as direct effects would spread support across cells.

**Bottom row — cancellation.** SNP 1 has effect 1 in both types, while SNP 2 has effect −1 in type 2. LD reduces the marginal effect of SNP 1 in type 2 to $1-0.8=0.2$. A shared joint effect can therefore look much weaker in one type.

This explains my concern about learning sharing patterns directly from marginal statistics. The comparison should specify whether a method accounts for LD **while learning the sharing prior**, during fine-mapping, or both. Conclusions about mvSuSiE should also specify the prior configuration being evaluated.

## 3. CASE model, following my handwritten derivation

Hold one gene fixed. The unit of expression measurement is a **donor within a cell type**; the target is a **SNP's effect across cell types**.

| Symbol | Object and dimension |
|---|---|
| $c$, $N_c$ | Cell type; number of donors observed in that type |
| $M$, $C$ | Number of cis SNPs; number of cell types |
| $y_c$, $X_c$ | Adjusted, standardized expression; donor-by-SNP genotype matrix |
| $B$ | $M\times C$ joint-effect matrix; row $B_j$ describes SNP $j$ across types |
| $\widehat B$ | $M\times C$ marginal association estimates |
| $R$ | $M\times M$ LD matrix: dependence across SNPs |
| $V$ | $C\times C$ sample-adjusted covariance: dependence across cell-type statistics |
| $U_t$, $\pi_t$ | Effect covariance for sharing pattern $t$; its mixture weight |

### Expression regression → summary-statistic likelihood

The underlying regressions are

$$
y_c=X_c\beta_c+\epsilon_c,\qquad B=(\beta_1,\ldots,\beta_C).
$$

Using the notation in my notes,

$$
\widehat B=\left(\frac{X_1^\top y_1}{N_1},\ldots,
\frac{X_C^\top y_C}{N_C}\right),\qquad
\widehat B\mid B\ \dot\sim\ \operatorname{MN}(RB,R,V).
$$

The dot indicates a summary-statistic approximation. A common LD matrix assumes sufficiently comparable genotype correlation structures across the cell-type samples. The matrix-normal form can be read as

$$
\operatorname{Cov}(\widehat B_{jc},\widehat B_{kd}\mid B)
\approx R_{jk}V_{cd}.
$$

Thus $R$ describes correlations across SNP rows, while $V$ describes correlations across cell-type columns. This distinction resolves two different sources of apparent sharing.

For overlapping donors, the sample-size adjustment in my notes is

$$
V_{cd}=V_{y,cd}\frac{N_{cd}}{N_cN_d}.
$$

Here $N_{cd}$ counts shared donors and $V_y$ is the phenotypic covariance used by the approximation. With standardized expression, $V_{cc}=1/N_c$. With no donor overlap, the off-diagonal sampling term is zero under the model assumptions. For example, $N_c=N_d=800$, $N_{cd}=600$, and $V_{y,cd}=0.4$ give $V_{cd}=0.000375$.

### Why a mixture of covariance matrices?

My notes write a mixture prior for one SNP's effect vector:

$$
B_j\sim\sum_{t=1}^{T}\pi_t N_C(0,U_t),\qquad
\sum_t\pi_t=1.
$$

Each matrix represents a possible pattern. For three cell types, examples are

$$
U_{\rm type\ 1}=s^2
\begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix},\qquad
U_{\rm types\ 1,2}=s^2
\begin{pmatrix}1&\rho&0\\\rho&1&0\\0&0&0\end{pmatrix}.
$$

The first allows an effect only in type 1. The second allows effects in types 1 and 2, with correlation $\rho$ and potentially different realized magnitudes. Here $|\rho|\le1$ ensures a valid covariance. An all-zero matrix supplies a point mass at the null effect vector.

A zero diagonal variance fixes that cell type's effect at zero under the pattern. Positive diagonal entries allow effects; off-diagonal entries describe how they co-vary. A zero off-diagonal alone means uncorrelated effects under that component. It does **not** imply that either cell type has zero effect. This is the distinction behind the “why $U_t$?” question in my notes.

### Fitting → posterior support

The model structure in my notes separates three operations:

1. Estimate sampling covariance $V$. For weak-signal SNPs $H$, the moment $|H|^{-1}\sum_{j\in H}\widehat B_j^\top\widehat B_j$ motivates a covariance estimate, followed by sample-size scaling. This relies on weak signals contributing little genetic mean.
2. Learn $\pi_t$ and $U_t$ using Monte Carlo EM: sample latent effects/patterns under the LD-aware likelihood, then update the prior parameters.
3. With the fitted prior, sample the posterior and estimate

$$
\operatorname{PIP}_{jc}=P(B_{jc}\ne0\mid\widehat B,R,\widehat V,\widehat\pi,\widehat U).
$$

The PIP targets a SNP–cell-type pair. A high value supports a modeled effect; its interpretation depends on the candidate variants, LD, prior, and sampling approximation. Very strong LD can leave multiple SNPs plausible even when a locus is clearly associated.

The [authors' package](https://github.com/leaffur/CASE) provides the implementation. I checked its [input interface](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE.R): summary statistics must align with LD; a sample-size vector assumes pairwise overlap equal to the smaller sample size. An explicit overlap matrix represents other designs. The optional input covariance defaults to independence, so its construction needs an explicit decision.

### Credible sets: what probability is being summarized?

The [inspected utility code](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE_uitility.R) groups variants using defaults of cumulative PIP at least 0.95 and pairwise absolute LD correlation at least 0.5. This clarifies the unsigned LD threshold in my notes.

There is a useful probability distinction. In a multi-causal model,

$$
\sum_{j\in A}\operatorname{PIP}_{jc}
=E[\text{number of active SNPs in }A\mid\text{data}],
$$

whereas set coverage concerns $P(\text{at least one active SNP in }A\mid\text{data})$. Equality holds when at most one SNP in the set can be active. I therefore treat the reported sets according to CASE's construction and assess their empirical coverage, rather than deriving a universal coverage guarantee from the PIP sum alone.

## 4. Offline investigation: what I checked

This entry combines my photographed notes, a source-code reading, and a small local calculation. **A fresh CASE fit on OneK1K was not run for this entry.**

| Check | Result and interpretation |
|---|---|
| Multiply $R B$ for the two sketches | Reproduces apparent sharing and cancellation in Figure 1 |
| Recover $B=R^{-1}\widehat B$ in the noise-free example | Recovers the specified effects; the matrix has condition number 9 |
| Inspect the package interface | Sample overlap and covariance inputs require attention |
| Inspect credible-set code | Confirms use of absolute LD correlation in the implementation |
| Compare with my September 28 report | Recovers the preliminary gene-overlap results below |

The toy inversion is an algebra check. With noisy data, nearly singular LD amplifies uncertainty; regularization and posterior inference become important. The heatmap illustrates the problem CASE addresses, without measuring CASE's performance.

[Download the reproducible Python calculation]({{base}}/assets/case-study/ld_example.py) · [Vector figure]({{base}}/assets/case-study/ld-effects.svg)

## 5. CIGMA–CASE overlap: pipeline and preliminary results

My September 28 progress report contains an existing downstream comparison. The workflow is:

**gene-level CIGMA calls → intersect tested gene IDs → join CASE credible-set calls → align broad lineages → compare detection and sharing**.

The local CIGMA analysis used 785 donors and seven finer cell types. CASE used eight broader lineages with differing donor availability. For the cell comparison, the report aligned **B, CD4, CD8, and NK**. The gene-level universe was restricted to 8,285 genes tested by both pipelines.

| Within the common tested universe | CASE eGene | CASE without an eGene call |
|---|---:|---:|
| CIGMA eGene | 316 | 0 |
| CIGMA without an eGene call | 3,524 | 4,445 |

The report records 3,840 CASE eGenes in that universe. The table is reconstructed from these reported totals: $3840-316=3524$, and $8285-3840=4445$. Among CIGMA non-calls, CASE identifies $3524/7969=44.2\%$ as eGenes.

This is **agreement between discovery lists**. Different preprocessing, donor eligibility, expression scales, and testing targets prevent interpreting the table as a controlled estimate of either method's precision or power. The original per-gene comparison tables remain necessary for a fresh audit; this entry transcribes the report and checks its arithmetic.

Among CASE eGenes in the respective CIGMA groups, the report gives credible sets in all four matched lineages for approximately **76% of cs-eGenes** and **97% of shared-only eGenes**. A gene can have supported effects in every lineage and still show substantial differences in effect magnitude. Finer T-cell or B-cell subtype differences are also collapsed by this alignment. Presence across broad lineages and homogeneity of effects answer different questions.

## 6. Downstream biological interpretation

Both approaches connect regulatory variation to biological annotations and disease genetics, using different starting objects.

| Starting object | Downstream comparison | Interpretation |
|---|---|---|
| CIGMA gene-level specificity | Gene properties and mapped genomic annotations | Which kinds of genes show heterogeneous genetic regulation? |
| CASE eGenes or variant sets | Cell markers, pathways, functional SNP annotations | Which biological contexts support the fine-mapped signals? |
| SNP annotations and GWAS summary statistics | Stratified LD-score regression | Is trait heritability concentrated in an annotation? |
| Matched eQTL and GWAS locus data | Additional colocalization analysis | Is a shared causal signal supported at this locus? |

CASE reported enhancer enrichment for cell-type-specific variants. Its eQTL annotations were heritability-enriched for five autoimmune diseases versus all SNPs and marginal eQTLs; differences versus SuSiE/mvSuSiE were statistically nonsignificant. Universally shared sets generally had the highest enrichment in its sharing-category comparison. [CASE Figure 6 and results](https://www.nature.com/articles/s41467-026-72176-3)

CIGMA reported enrichment of cell-type-specific regulation in complex-trait heritability. [CIGMA paper](https://www.nature.com/articles/s41586-026-10577-6) These observations use different annotations, selection rules, and estimands. Comparing them requires matching the SNP universe, gene-to-SNP mapping, expression scale, and power-related filters.

For this diary, the local overlap supports a precise interpretation: **a heterogeneous gene-level genetic signal can coexist with fine-mapping support across several lineages**. The next level of evidence is the aligned effect profile and its uncertainty for each locus.

---

*Study record: September 30, 2026. Equations organized from my handwritten CASE notes; primary references are Chen & Dahl's CIGMA paper and Lin et al.'s CASE paper. Local overlap numbers come from my September 28 progress report. The LD figure is an original, executed teaching calculation. Package interface and utility code inspected at commit `13f4fc8`.*
