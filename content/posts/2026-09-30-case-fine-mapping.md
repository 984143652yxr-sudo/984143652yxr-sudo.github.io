---
{
  "slug": "case-fine-mapping",
  "title": "CIGMA and CASE: variance components, fine-mapping, and cellular context",
  "date": "2026-09-30",
  "date_label": "September 30, 2026",
  "status": "Statistical study note",
  "reading_time": "13 min read",
  "summary": "Gaussian cis-SNP effects → joint effect matrices → credible sets and eGenes → overlap and biological enrichment."
}
---

## 1. CIGMA: genetic variance across cell types

Fix one gene and retain $M$ SNPs in its cis region. Let $Y_{ic}$ be expression for donor $i$ in cell type $c$, and let $Z_{ij}$ be standardized genotype at SNP $j$. A simplified shared-plus-specific model is

$$
Y_{ic}=\mu_c+W_i^\top\gamma_c+
\underbrace{\sum_{j=1}^{M}Z_{ij}(\alpha_j+\eta_{jc})}_{a_{ic}:\ \text{genetic contribution}}+r_{ic},
$$

$$
\alpha_j\sim N\left(0,\frac{\sigma^2_{\rm shared}}{M}\right),\qquad
\eta_{jc}\sim N\left(0,\frac{v_c}{M}\right).
$$

$\alpha_j$ is shared across cell types; $\eta_{jc}$ is a cell-type deviation. The independent Gaussian coefficients extend over **all retained cis SNPs**. The $1/M$ scaling makes the parameters describe aggregate genetic variance. With $K=ZZ^\top/M$,

$$
\operatorname{Cov}(a_{ic},a_{\ell d}\mid Z)
=K_{i\ell}\left[\sigma^2_{\rm shared}+\mathbf1(c=d)v_c\right].
$$

| Parameter or test | Target |
|---|---|
| $\sigma^2_{\rm shared}$ | Shared genetic variance |
| $v_c$ | Additional genetic variance specific to type $c$ in this construction |
| $H_0:v_1=\cdots=v_C=0$ | A gene-level test of cell-type specificity |
| $\bar v/(\sigma^2_{\rm shared}+\bar v)$ | Aggregate specificity, where $\bar v=C^{-1}\sum_c v_c$ |

The fitted CIGMA specification also accounts for residual dependence and finite-cell measurement uncertainty. The equations above isolate its genetic component. [CIGMA: Chen & Dahl (2026)](https://www.nature.com/articles/s41586-026-10577-6)

### From a variance parameter to a causal-variant question

Integrating over the Gaussian coefficients produces a likelihood for variance components. A nondegenerate Gaussian prior assigns probability zero to an exactly zero coefficient; it supplies no SNP-level inclusion indicator. CIGMA therefore characterizes **how much genetic regulation differs across cell types**, while its standard outputs leave individual causal-SNP localization unresolved.

A gene can have $v_c>0$ even when it has effects in every cell type: differences in magnitude are sufficient. Several SNP-effect configurations can also produce similar aggregate covariance. Fine-mapping introduces a sparse, variant-level model to distinguish plausible configurations using LD.

**Quantitative motivation.** Published CIGMA results estimated cis specificity at 31.1%, compared with 59.4% for trans regulation; a bulk-like expression aggregate gave about 5% cis specificity. Its OneK1K analysis detected 193 cs-eGenes among 10,288 genes. These figures quantify variance and gene-level heterogeneity, respectively. [CIGMA results](https://www.nature.com/articles/s41586-026-10577-6)

## 2. CASE: one effect matrix for all cell types

For the same gene, organize all cell types jointly:

$$
y_c=X_c\beta_c+\epsilon_c,\qquad c=1,\ldots,C,
$$

$$
B=(\beta_1,\ldots,\beta_C)
=\begin{pmatrix}
\beta_{11}&\cdots&\beta_{1C}\\
\vdots&\ddots&\vdots\\
\beta_{M1}&\cdots&\beta_{MC}
\end{pmatrix}\in\mathbb R^{M\times C}.
$$

**Rows are SNPs; columns are cell types.** Column $\beta_c$ contains the cis effects for type $c$. Row $B_j=(\beta_{j1},\ldots,\beta_{jC})$ contains one SNP's effects across all types. Different types can have different numbers of donors, so the $X_c$ matrices can have different row counts.

| Symbol | Meaning |
|---|---|
| $y_c$, $X_c$ | Adjusted standardized expression; standardized donor-by-SNP genotypes |
| $N_c$, $N_{cd}$ | Donors in type $c$; donors overlapping types $c,d$ |
| $B$, $\widehat B$ | Joint effects; marginal association estimates, both $M\times C$ |
| $R$ | $M\times M$ LD matrix across SNPs |
| $V$ | $C\times C$ sampling covariance across cell-type statistics |
| $U_t$, $\pi_t$ | Effect covariance for sharing pattern $t$; pattern probability |

### LD separates marginal associations from joint effects

For standardized variables,

$$
\widehat b_c=X_c^\top y_c/N_c,\qquad
E[\widehat b_c\mid\beta_c]\approx R\beta_c.
$$

<figure class="study-figure" style="max-width:600px">
<a href="{{base}}/assets/case-study/ld-simple.svg"><img src="{{base}}/assets/case-study/ld-simple.svg" alt="Two SNPs each affect a different cell type, but LD creates marginal associations in both types"></a>
<figcaption>Figure 1. A two-SNP example with LD correlation 0.8. Joint effects are sparse across types; marginal associations spread across both columns. Values are illustrative.</figcaption>
</figure>

The first SNP's marginal effect in type 2 is $0+0.8\times1=0.8$, although its joint effect is zero. Opposing joint effects can instead cancel. The statistical task is to infer $B$ from LD-mixed observations $\widehat B$.

### Summary-statistic likelihood

CASE uses the approximate matrix-normal likelihood

$$
\widehat B\mid B\ \dot\sim\ \operatorname{MN}(RB,R,V),
\qquad
\operatorname{Cov}(\widehat B_{jc},\widehat B_{kd}\mid B)
\approx R_{jk}V_{cd}.
$$

$R$ models dependence across SNPs, and $V$ models dependence across cell-type statistics. The sample-overlap adjustment is

$$
V_{cd}=V_{y,cd}\frac{N_{cd}}{N_cN_d},\qquad V_{cc}=1/N_c
\quad\text{for standardized expression}.
$$

$V_y$ is the phenotypic covariance used by the approximation. A common $R$ requires comparable LD across the sampled donor groups. [CASE model: Lin et al. (2026)](https://www.nature.com/articles/s41467-026-72176-3)

### A sparse prior on each SNP row

Introduce a latent sharing pattern $z_j$:

$$
z_j\sim\operatorname{Categorical}(\pi_0,\ldots,\pi_T),\qquad
B_j\mid z_j=t\sim N_C(0,U_t),\qquad U_0=0.
$$

The null component sets the entire row to zero. Other components allow effects in selected cell types. For three types, two illustrative patterns are

$$
U_{\{1\}}=s^2\begin{pmatrix}1&0&0\\0&0&0\\0&0&0\end{pmatrix},\qquad
U_{\{1,2\}}=s^2\begin{pmatrix}1&\rho&0\\\rho&1&0\\0&0&0\end{pmatrix},
\quad |\rho|\le1.
$$

A zero diagonal fixes that cell's effect at zero. Positive diagonals allow nonzero effects; off-diagonals describe their correlation. Shared activity can have different effect magnitudes. General fitted matrices can also have unequal variances across active types.

The inference sequence is

$$
(\widehat B,R,N)\ \longrightarrow\ \widehat V
\ \longrightarrow\ (\widehat\pi_t,\widehat U_t)
\ \longrightarrow\ P(B\mid\text{data})
\ \longrightarrow\operatorname{PIP}_{jc}.
$$

Monte Carlo EM learns the sharing mixture under the LD-aware likelihood. Posterior sampling then estimates

$$
\operatorname{PIP}_{jc}=P(\beta_{jc}\ne0\mid\widehat B,R,\widehat V,\widehat\pi,\widehat U).
$$

The model borrows information across cell types through $U_t$, while the likelihood accounts for LD through $R$. Implementation details and input conventions are given in the [CASE source](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE.R).

## 3. From PIPs to credible sets, eQTLs, and eGenes

### Cell-type-wise credible-set construction

For gene $g$ and type $c$, let $A$ be a candidate SNP set. CASE's default set criteria are

$$
\sum_{j\in A}\operatorname{PIP}_{jc}\ge0.95,\qquad
\min_{j\ne k\in A}|R_{jk}|\ge0.5.
$$

The [implementation](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE_models.R) scans SNPs in descending PIP order, reports singleton sets for PIP at least 0.95, and otherwise searches LD-compatible candidates for sufficient cumulative PIP. Selected variants are flagged before continuing. The [search utility](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE_uitility.R) returns the first valid set under its ordering; this is an operational search rather than a demonstrated global minimum-cardinality solution.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/credible-set-example.svg"><img src="{{base}}/assets/case-study/credible-set-example.svg" alt="A worked example in which two correlated SNPs accumulate PIP 0.96 in cell type A, while type B has insufficient posterior support"></a>
<figcaption>Figure 2. Illustrative PIPs for one gene. In type A, SNPs 1 and 2 form a reported set; type B has insufficient support. These numbers demonstrate the calling rule.</figcaption>
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

## 4. CIGMA–CASE overlap

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

## 5. Potential limitations of each model

| Issue | CIGMA | CASE |
|---|---|---|
| Inferential target | Estimates aggregate genetic variance; causal SNP identities remain unresolved | Prioritizes SNP–cell-type effects; strong LD can leave several variants plausible |
| Effect architecture | Gaussian variance model averages over the retained SNP set | Sparse mixture depends on the available sharing patterns and their estimation |
| Weak signals | Gene-level variance estimates can be imprecise | Low PIP or a missing set may reflect limited power, especially in rare types |
| Cellular resolution | Pseudobulk averages over within-type state variation | Cell-type summaries also average over within-type states |
| Dependence and measurement | Requires appropriate residual and sampling covariance | Requires aligned LD, sample-overlap information, and sampling covariance |
| Scientific interpretation | Specific variance supports heterogeneous genetic regulation | PIPs and sets provide model-based causal candidates; functional validation remains a separate step |

For both approaches, sample size, expression scale, and annotation choices affect interpretation. A comparison of significant lists alone combines these influences with differences between the models.

## 6. Downstream biological interpretation: construction and models

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/downstream-process.svg"><img src="{{base}}/assets/case-study/downstream-process.svg" alt="Three downstream routes: gene-set enrichment, SNP-annotation heritability enrichment, and locus-level colocalization"></a>
<figcaption>Figure 3. Three analysis routes from regulatory results to biological evidence. Colocalization is a possible locus-level follow-up; it requires an additional model.</figcaption>
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

**Data and computation.** The gene-overlap counts are transcribed from the September 28 OneK1K progress report and have not been re-estimated here. The LD and credible-set diagrams are teaching examples. [Figure-generation code]({{base}}/assets/case-study/make_diagrams.py). CASE source inspection refers to commit `13f4fc8`; no new OneK1K CASE fit is reported.
