---
{
  "slug": "case-fine-mapping",
  "title": "CIGMA and CASE: variance components, fine-mapping, and cellular context",
  "date": "2026-09-30",
  "date_label": "September 30, 2026",
  "status": "Statistical study note",
  "reading_time": "17 min read",
  "summary": "One donor-by-cell-type regression, two distributions for SNP effects: covariance estimation with CIGMA and posterior fine-mapping with CASE."
}
---

## 1. One regression framework for a fixed gene

Both methods start from **donor-by-cell-type expression and cis-SNP genotypes for one gene**. Cells are aggregated within donor and cell type; the genetic observation is a donor. Fix the gene and suppress its index throughout.

After the chosen expression transformation and appropriate covariate adjustment, write

$$
y_{ic}=\sum_{j=1}^{M}x_{ij}\beta_{jc}+e_{ic}.
$$

For a common set of $N$ donors measured in $C$ cell types,

$$
\boxed{Y=XB+E},\qquad
Y\in\mathbb R^{N\times C},\quad
X\in\mathbb R^{N\times M},\quad
B\in\mathbb R^{M\times C}.
$$

$$
B=(\beta_1,\ldots,\beta_C)
=\begin{pmatrix}
\beta_{11}&\cdots&\beta_{1C}\\
\vdots&\ddots&\vdots\\
\beta_{M1}&\cdots&\beta_{MC}
\end{pmatrix},\qquad
b_j=(\beta_{j1},\ldots,\beta_{jC})^\top.
$$

**A column of $B$ contains all cis-SNP effects in one cell type. A row contains one SNP's effects across all cell types.** The column vector $b_j$ is that row transposed. Both CIGMA and CASE specify a distribution for these same SNP-effect vectors.

| Symbol | Meaning | Dimension |
|---|---|---|
| $Y$, $y_c$ | Donor-level expression across types; column for type $c$ | $N\times C$; $N\times1$ |
| $X$ | Centered, standardized genotypes for retained cis SNPs | $N\times M$ |
| $B$, $b_j$ | Joint SNP effects; one SNP's effect vector across types | $M\times C$; $C\times1$ |
| $E$ | Remaining biological and measurement variation | $N\times C$ |
| $K=XX^\top/M$ | Genetic similarity **between donors** | $N\times N$ |
| $R=X^\top X/N$ | LD correlation **between SNPs** | $M\times M$ |
| $\Omega$ | Aggregate genetic covariance across cell types | $C\times C$ |
| $\widehat B$, $V$ | Marginal SNP associations; their across-type sampling covariance | $M\times C$; $C\times C$ |

$K$ and $R$ summarize the same genotype matrix along different axes. This explains why a variance-component method uses donor relatedness while fine-mapping uses SNP LD.

With different donor sets across types, use $y_c=X_c\beta_c+e_c$ and retain the observed donor identities. The balanced matrix form above is a common notation, rather than a requirement to discard incomplete donors. Expression normalization, covariate treatment, and residual specifications must still follow each method; identical symbols do not make their preprocessing identical.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/unified-framework.svg"><img src="{{base}}/assets/case-study/unified-framework.svg" alt="A common donor-by-cell-type regression branches into CIGMA covariance estimation and CASE sparse posterior fine-mapping"></a>
<figcaption>One response and effect matrix, two inferential targets. The Gaussian branch estimates aggregate covariance; the sparse-mixture branch retains SNP-level inclusion uncertainty.</figcaption>
</figure>

## 2. The difference is the distribution of the same SNP-effect vector

### CIGMA: one structured Gaussian distribution

In the shared-plus-specific construction,

$$
\beta_{jc}=\alpha_j+\eta_{jc},\qquad
\alpha_j\sim N\left(0,\frac{\sigma^2_{\rm shared}}{M}\right),\qquad
\eta_{jc}\sim N\left(0,\frac{v_c}{M}\right).
$$

Take the shared and specific coefficients independent across SNPs, and the specific coefficients independent across types. Then the entire SNP row has distribution

$$
\boxed{b_j\sim N_C(0,\Omega/M)},\qquad
\Omega=\sigma^2_{\rm shared}\mathbf1\mathbf1^\top+
\operatorname{diag}(v_1,\ldots,v_C).
$$

For three types,

$$
\Omega=
\begin{pmatrix}
\sigma^2_{\rm shared}+v_1&\sigma^2_{\rm shared}&\sigma^2_{\rm shared}\\
\sigma^2_{\rm shared}&\sigma^2_{\rm shared}+v_2&\sigma^2_{\rm shared}\\
\sigma^2_{\rm shared}&\sigma^2_{\rm shared}&\sigma^2_{\rm shared}+v_3
\end{pmatrix}.
$$

All retained cis SNPs enter the working random-effect distribution. Its parameters summarize genetic covariance. The Gaussian construction provides no SNP-specific activity indicator; this avoids requiring a list of individually detected variants before estimating aggregate variation. CIGMA also supports a more general cross-type genetic covariance than the shared-plus-diagonal form shown here. [CIGMA, Chen et al. (2026)](https://www.nature.com/articles/s41586-026-10577-6)

### CASE: a mixture of sparse sharing patterns

Keep the same $b_j$, but introduce a pattern indicator $z_j$:

$$
z_j\sim\operatorname{Categorical}(\pi_0,\ldots,\pi_T),\qquad
\boxed{b_j\mid z_j=t\sim N_C(0,U_t)},\qquad U_0=0.
$$

Each $U_t$ is a $C\times C$ covariance. Its zero rows and columns specify inactive cell types. Examples for three types are

$$
U_{\{1\}}=s^2\begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix},\qquad
U_{\{1,2\}}=s^2\begin{pmatrix}1&\rho&0\\\rho&1&0\\0&0&0\end{pmatrix},
\quad |\rho|\le1.
$$

The first permits an effect only in type 1; the second permits effects in types 1 and 2. The all-zero component fixes the entire SNP row to zero. Within an active subset, effects may have different magnitudes. CASE learns mixture weights and effect covariances and uses the fitted mixture for posterior fine-mapping. [CASE, Lin et al. (2026)](https://www.nature.com/articles/s41467-026-72176-3)

### Their exact second-moment connection

Assume independent, zero-mean SNP-effect vectors conditional on the prior parameters. Under CASE,

$$
\operatorname{Var}(b_j)=\sum_{t=0}^{T}\pi_tU_t,
\qquad
\Omega_{\rm CASE}=M\sum_{t=0}^{T}\pi_tU_t.
$$

Writing $A=XB$ for the genetic contribution, both constructions imply

$$
\boxed{\operatorname{Cov}(A_{ic},A_{\ell d}\mid X)
=K_{i\ell}\Omega_{cd}}.
$$

For CIGMA, $\Omega$ has its specified variance-component structure. For CASE, the same identity uses $\Omega_{\rm CASE}$. The covariance can agree even when the underlying distributions differ.

For example, compare a Gaussian prior with a sparse mixture:

$$
b_j\sim N_C(0,\Omega/M)
\quad\text{and}\quad
b_j\sim0.9\,\delta_0+0.1\,N_C(0,10\Omega/M).
$$

Both have covariance $\Omega/M$. The second assigns 90% prior probability to an entirely inactive SNP. The first uses a continuous effect distribution. Aggregate covariance alone cannot distinguish these two configurations. Fine-mapping adds distributional assumptions and SNP-level evidence to infer activity.

Here $\delta_0$ denotes a point mass at the zero vector. This is an algebraic comparison of priors, not a claim that the complete algorithms are interchangeable.

## 3. Different inference targets lead to different tasks

### CIGMA: integrate SNP effects and estimate covariance components

Stack expression by cell type. With Gaussian residuals independent of the SNP effects,

$$
\operatorname{vec}(Y)\sim
N\left(0,\Omega\otimes K+\Sigma_E\right),
$$

where $\Sigma_E$ represents the method's residual and measurement covariance. For example,

$$
\operatorname{Cov}(Y_{ic},Y_{\ell d}\mid X)
=K_{i\ell}\Omega_{cd}+\operatorname{Cov}(e_{ic},e_{\ell d}).
$$

CIGMA's HE approach estimates variance components from second moments, with jackknife inference. A schematic moment-fitting objective, after appropriate fixed-effect projection, is

$$
\widehat\theta=\arg\min_\theta
\left\|yy^\top-\Sigma_0-\sum_q\theta_qA_q\right\|_F^2.
$$

Here $y$ is stacked adjusted expression, $A_q$ are known covariance-component matrices, and $\Sigma_0$ collects specified covariance offsets. Genetic and residual components enter the fit together. This equation explains the estimation principle; the implementation supplies its projection and covariance details. [CIGMA methods](https://www.nature.com/articles/s41586-026-10577-6)

Write $\bar v=C^{-1}\sum_c v_c$ for the mean cell-type-specific variance. The outputs concern the **gene's genetic architecture**:

$$
\widehat\sigma^2_{\rm shared},\quad\widehat v_c,\quad
\frac{\widehat{\bar v}}{\widehat\sigma^2_{\rm shared}+\widehat{\bar v}},
\qquad H_0:v_1=\cdots=v_C=0.
$$

Gaussian models can also yield shrunken effect predictions if that calculation is added. A continuous prior alone does not supply causal-inclusion probabilities or a sparse causal set. CIGMA's standard target is variance estimation and testing.

### CASE: retain SNP effects and infer their posterior distribution

Marginal association estimates satisfy

$$
\widehat\beta_c=X_c^\top y_c/N_c,\qquad
E[\widehat\beta_c\mid\beta_c]\approx R\beta_c.
$$

<figure class="study-figure" style="max-width:600px">
<a href="{{base}}/assets/case-study/ld-simple.svg"><img src="{{base}}/assets/case-study/ld-simple.svg" alt="Two SNPs each affect a different cell type, but LD creates marginal associations in both types"></a>
<figcaption>LD mixes the rows of B. Each SNP acts in one type in this example, but has marginal associations in both. Values are illustrative.</figcaption>
</figure>

With $\widehat B=(\widehat\beta_1,\ldots,\widehat\beta_C)$,

$$
\widehat B\mid B\ \dot\sim\ \operatorname{MN}(RB,R,V),\qquad
\operatorname{Cov}(\widehat B_{jc},\widehat B_{kd}\mid B)
\approx R_{jk}V_{cd}.
$$

$V$ is the sampling covariance of the marginal statistics. It is distinct from the biological effect covariance $\Omega$. CASE's approximation uses

$$
V_{cd}=V_{y,cd}\frac{N_{cd}}{N_cN_d},\qquad V_{cc}=1/N_c
\quad\text{for standardized expression}.
$$

The donor overlap $N_{cd}$ and phenotypic covariance $V_y$ determine the adjustment. The common-$R$ approximation also requires comparable LD across donor groups.

Let $\psi=\{\pi_t,U_t\}$. Monte Carlo EM learns $\psi$ by alternating posterior sampling of $(B,z)$ with updates of the prior parameters. With the fitted prior,

$$
p(B,z\mid\widehat B,R,\widehat V,\widehat\psi)
\propto p(\widehat B\mid B,R,\widehat V)
\prod_{j=1}^{M}p(b_j,z_j\mid\widehat\psi).
$$

For $S$ retained posterior draws,

$$
\operatorname{PIP}_{jc}=P(\beta_{jc}\ne0\mid\text{data}),\qquad
\widehat{\operatorname{PIP}}_{jc}
\approx\frac1S\sum_{s=1}^{S}\mathbf1(\beta_{jc}^{(s)}\ne0).
$$

LD couples the SNP rows in the likelihood; the sharing mixture couples cell types within each row. The output is **variant-by-cell-type support**, followed by credible sets. [CASE source and input conventions](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE.R)

### Match the output to the scientific question

| Scientific task | Target in the common notation | Appropriate output |
|---|---|---|
| Quantify shared and heterogeneous cis regulation | Structure and magnitude of $\Omega$ | CIGMA variance estimates and uncertainty |
| Test whether a gene's genetic effects differ across types | $v_1=\cdots=v_C=0$ in the shared-plus-specific model | CIGMA gene-level specificity test |
| Prioritize an experimentally testable variant | Which rows/entries of $B$ are active? | CASE PIPs and credible sets |
| Identify the cell types supporting a variant | Pattern of nonzero entries in $b_j$ | CASE SNP-level posterior support |
| Estimate whether regulation is detectable for a gene in a type | At least one supported signal in column $\beta_c$ | CASE eGene call |
| Relate regulatory architecture to disease heritability | Gene or SNP annotations built from these outputs | Annotation-specific downstream analysis |

**Shared activity and equal effects are different targets.** The vector $(0.1,0.3,0.7)$ is active in every type and heterogeneous in magnitude. CASE can support its activity pattern while CIGMA can detect aggregate heterogeneity. Neither a CASE non-call nor an imprecise variance estimate establishes absence of regulation.

## 4. From PIPs to credible sets, eQTLs, and eGenes

### Cell-type-wise credible-set construction

For gene $g$ and type $c$, let $A$ be a candidate SNP set. CASE's default set criteria are

$$
\sum_{j\in A}\operatorname{PIP}_{jc}\ge0.95,\qquad
\min_{j\ne k\in A}|R_{jk}|\ge0.5.
$$

The [implementation](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE_models.R) scans SNPs in descending PIP order, reports singleton sets for PIP at least 0.95, and otherwise searches LD-compatible candidates for sufficient cumulative PIP. Selected variants are flagged before continuing. The [search utility](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE_uitility.R) returns the first valid set under its ordering; this is an operational search rather than a demonstrated global minimum-cardinality solution.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/credible-set-example.svg"><img src="{{base}}/assets/case-study/credible-set-example.svg" alt="A worked example in which two correlated SNPs accumulate PIP 0.96 in cell type A, while type B has insufficient posterior support"></a>
<figcaption>Illustrative PIPs for one gene. In type A, SNPs 1 and 2 form a reported set; type B has insufficient support. These numbers demonstrate the calling rule.</figcaption>
</figure>

The output for one gene is a collection of sets **for each cell type**:

$$
\mathcal C_{gc}=\{A_{gc1},\ldots,A_{gcL_{gc}}\}.
$$

Several sets may represent several signals. SNPs within a set remain competing or jointly plausible candidates. The term “putative causal eQTLs” refers to variants prioritized in these sets; set membership alone does not establish every member as causal.

### The eGene calling rule

A gene has an eGene call in cell type $c$ when at least one credible set is reported:

$$
E_{gc}=\mathbf1(L_{gc}\ge1),\qquad
E_g=\mathbf1\left(\sum_{c=1}^{C}E_{gc}\ge1\right).
$$

The vector $(E_{g1},\ldots,E_{gC})$ describes the **gene's detected cell-type pattern**. Different cell types can support different variants for the same gene. Gene-level sharing therefore does not by itself establish sharing of the same causal SNP.

### Interpreting “95%”

The posterior probability of a set containing an active variant is

$$
P\left(\sum_{j\in A}\mathbf1(\beta_{jc}\ne0)\ge1\mid\text{data}\right).
$$

The PIP sum instead equals the posterior expected number of active variants. They coincide under an at-most-one-active-variant assumption within the set. In general multi-causal settings, CASE's nominal set rule should be assessed through coverage calibration.

**Published result.** CASE identified 5,057 eGenes among 11,704 candidate genes in OneK1K, 13.0% more than mvSuSiE and 32.9% more than SuSiE in that comparison. Its simulations varied sample size, heritability, and sharing patterns, including null cell types to test over-sharing. [CASE results](https://www.nature.com/articles/s41467-026-72176-3)

## 5. CIGMA–CASE overlap

A gene-level comparison requires a common tested universe:

$$
\mathcal G=\mathcal G_{\rm tested,CIGMA}\cap\mathcal G_{\rm tested,CASE}
\quad\longrightarrow\quad
\text{join gene calls}\quad\longrightarrow\quad
\text{compare cell-type patterns}.
$$

The preliminary OneK1K comparison used 8,285 common genes. For lineage-level interpretation, the seven finer CIGMA types were aligned to **B, CD4, CD8, and NK** in CASE. The CIGMA eGene list combines its shared and specific discovery categories.

| Gene calls within the common universe | CASE eGene | CASE non-call |
|---|---:|---:|
| CIGMA eGene | 316 | 0 |
| CIGMA non-call | 3,524 | 4,445 |

Thus CASE recovered all 316 CIGMA gene calls and identified an additional 3,524 genes. Among CASE eGenes within the CIGMA categories, approximately **76% of cs-eGenes** and **97% of shared-only eGenes** had credible sets in all four matched lineages.

The overlap concerns detection. CIGMA specificity concerns effect heterogeneity. For example,

$$
(\beta_{j1},\beta_{j2},\beta_{j3})=(0.1,0.3,0.7)
$$

is compatible with activity in all three types and substantial differences in magnitude. Donor eligibility, cell grouping, expression normalization, and significance criteria differ between the compared pipelines, so these counts describe concordance rather than a controlled power benchmark.

## 6. Potential limitations of each model

| Issue | CIGMA | CASE |
|---|---|---|
| Inferential target | Estimates aggregate genetic variance; causal SNP identities remain unresolved | Prioritizes SNP–cell-type effects; strong LD can leave several variants plausible |
| Effect architecture | Gaussian variance model averages over the retained SNP set | Sparse mixture depends on the available sharing patterns and their estimation |
| Weak signals | Gene-level variance estimates can be imprecise | Low PIP or a missing set may reflect limited power, especially in rare types |
| Cellular resolution | Pseudobulk averages over within-type state variation | Cell-type summaries also average over within-type states |
| Dependence and measurement | Requires appropriate residual and sampling covariance | Requires aligned LD, sample-overlap information, and sampling covariance |
| Scientific interpretation | Specific variance supports heterogeneous genetic regulation | PIPs and sets provide model-based causal candidates; functional validation remains a separate step |

For both approaches, sample size, expression scale, and annotation choices affect interpretation. A comparison of significant lists alone combines these influences with differences between the models.

## 7. Downstream biological interpretation: construction and models

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/downstream-process.svg"><img src="{{base}}/assets/case-study/downstream-process.svg" alt="Three downstream routes: gene-set enrichment, SNP-annotation heritability enrichment, and locus-level colocalization"></a>
<figcaption>Three analysis routes from regulatory results to biological evidence. Colocalization is a possible locus-level follow-up; it requires an additional model.</figcaption>
</figure>

### A. Gene calls → marker or pathway enrichment

Define a tested-gene universe $\mathcal G$, a selected set $D$ (for example, cell-type-specific eGenes), and a marker/pathway set $F$. Construct

$$
\begin{array}{c|cc}
& F & \mathcal G\setminus F\\\hline
D & a & b\\
\mathcal G\setminus D & c & d
\end{array}
\qquad
\operatorname{OR}=\frac{ad}{bc}.
$$

Fisher's exact test assesses enrichment conditional on the margins; adjust across tested pathways or cell types. The universe should contain eligible tested genes. Matching on expression, cell abundance, or other detection-related properties can assess sensitivity to selection effects.

### B. SNP sets → annotations → trait-heritability enrichment

For a binary SNP annotation $A$, let $a_A(k)=1$ if SNP $k$ belongs to the selected category. For CASE, one possible category is the union of credible-set variants in selected cell types. For CIGMA, a declared gene-to-SNP mapping is needed to translate gene-level outputs into genomic annotations.

Construct annotation-specific LD scores and fit the stratified regression

$$
\ell_A(j)=\sum_k r_{jk}^2a_A(k),\qquad
E[\chi_j^2]\approx b+N_{\rm GWAS}\sum_A\tau_A\ell_A(j).
$$

Here $b$ is the intercept, $\tau_A$ is the annotation's conditional contribution, and the fitted model includes relevant baseline annotations. Report

$$
\operatorname{Enrichment}(A)=
\frac{h_A^2/h^2}{M_A/M_{\rm all}}.
$$

This compares the annotation's heritability share with its SNP share. Overlapping annotations require joint modeling; $\tau_A$ and enrichment are different summaries. Use ancestry-compatible LD and uncertainty estimates. [S-LDSC documentation](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability)

**Published findings.** CASE reported enhancer enrichment for cell-type-specific variants and autoimmune-trait heritability enrichment relative to all SNPs and marginal eQTLs. Differences versus SuSiE/mvSuSiE were statistically nonsignificant; universally shared categories generally showed the highest enrichment in its sharing comparison. [CASE Figure 6](https://www.nature.com/articles/s41467-026-72176-3) CIGMA reported enrichment of cell-type-specific regulation in complex-trait heritability. [CIGMA](https://www.nature.com/articles/s41586-026-10577-6) The annotation definitions and selection rules differ, motivating a harmonized comparison.

### C. eQTL locus + GWAS locus → colocalization

After allele harmonization and compatible fine-mapping, compare hypotheses

$$
H_{\rm shared}:\text{a shared causal variant},\qquad
H_{\rm distinct}:\text{different causal variants}.
$$

In a simplified one-causal-variant-per-trait calculation, with variant Bayes factors $\mathrm{BF}^{E}_j$ and $\mathrm{BF}^{G}_j$ and uniform location weights, the relative evidence contains

$$
\mathcal L_{\rm shared}\propto\sum_j\mathrm{BF}^{E}_j\mathrm{BF}^{G}_j,
\qquad
\mathcal L_{\rm distinct}\propto\sum_{j\ne k}\mathrm{BF}^{E}_j\mathrm{BF}^{G}_k.
$$

Posterior probabilities also require hypothesis priors, location-weight normalization, and the remaining association hypotheses. Multiple signals require a suitable multi-signal analysis. An intersection of credible sets can prioritize a locus for this calculation; colocalization quantifies the evidence for sharing, while mediation requires further assumptions. [Coloc model and hypotheses](https://chr1swallace.github.io/coloc/articles/a03_enumeration.html)

---

**Revision: October 1, 2026. Data and computation.** The gene-overlap counts are transcribed from the September 28 OneK1K progress report and have not been re-estimated here. The LD and credible-set diagrams are teaching examples. [Figure-generation code]({{base}}/assets/case-study/make_diagrams.py). CASE source inspection refers to commit `13f4fc8`; no new OneK1K CASE fit is reported.
