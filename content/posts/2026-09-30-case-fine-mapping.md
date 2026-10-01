---
{
  "slug": "case-fine-mapping",
  "title": "CIGMA and CASE: variance components, fine-mapping, and cellular context",
  "date": "2026-09-30",
  "date_label": "September 30, 2026",
  "status": "Statistical study note",
  "reading_time": "20 min read",
  "summary": "Shared genetic-effect notation, followed by the full CIGMA variance model and CASE fine-mapping model: residuals, priors, inference, and scientific targets."
}
---

## 1. Shared notation: one gene, donors, cell types, and SNP effects

Fix one gene throughout. Each expression observation summarizes cells from **one donor and one cell type**. Genotype varies between donors; a donor's inherited genotype is shared across their cell types.

| Symbol | Meaning |
|---|---|
| $i$, $c$, $j$ | Donor, cell type, and SNP indices |
| $C$; $N$; $N_c$ | Number of cell types; total donors; donors observed in type $c$ |
| $y_{ic}$ | Expression summary for this gene in donor $i$, type $c$, on the method's analysis scale |
| $x_{ij}$ | Standardized genotype of donor $i$ at SNP $j$ |
| $\beta_{jc}$ | Joint additive effect of SNP $j$ on expression in type $c$ |
| $\beta_c$ | All included SNP effects in **one cell type**: a column vector |
| $b_j$ | One SNP's effects across **all cell types**: a column vector of length $C$ |
| $B$ | Effect matrix: **SNP rows × cell-type columns**, with $B_{jc}=\beta_{jc}$ |

The two views of the effect matrix are

$$
B=(\beta_1,\ldots,\beta_C),\qquad
b_j=(\beta_{j1},\ldots,\beta_{jC})^\top.
$$

Thus $\beta_c$ follows a column down the SNPs, while $b_j$ follows a row across cell types. A distribution for $b_j$ is a joint distribution for several effects, one per cell type.

The notation identifies the biological quantities. Each method below specifies its own expression scale, variant set, residual model, and likelihood. Use $L$ for the number of SNPs in a CIGMA component and $M$ for the CASE cis region.

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/case-study/unified-framework.svg"><img src="{{base}}/assets/case-study/unified-framework.svg" alt="Shared donor, cell-type, and SNP notation leads to a CIGMA pseudobulk covariance model and a CASE summary-statistic fine-mapping model"></a>
<figcaption>The shared objects are donors, cell types, and genetic effects. The observation models and inferential targets are specified separately.</figcaption>
</figure>

## 2. CIGMA: partition genetic and residual variation

### Observation model

For one genotype component containing $L$ SNPs, the shared-plus-specific model can be written

$$
\boxed{
y_{ic}=m_{ic}
+\sum_{j=1}^{L}x_{ij}(\alpha_j+\eta_{jc})
+u_i+r_{ic}+\epsilon_{ic}.}
$$

Here $m_{ic}$ contains the cell-type mean and fitted covariate effects. The random terms have different roles:

| Term | Role | Variance |
|---|---|---|
| $\alpha_j$ | SNP effect shared across types | $\sigma_g^2/L$ |
| $\eta_{jc}$ | SNP's deviation in type $c$ | $v_c/L$ |
| $u_i$ | Residual donor effect shared across types | $\sigma_e^2$ |
| $r_{ic}$ | Residual donor effect specific to type $c$ | $w_c$ |
| $\epsilon_{ic}$ | Sampling error in the cell-type mean | $\delta_{ic}$ |

The first four are zero-mean Gaussian random effects in this construction. They are mutually independent, with independence across SNPs or donors as appropriate; the same $\alpha_j$ and $u_i$ recur across types. Residual donor variation includes unmodelled genetic and nongenetic contributions.

For $n_{ic}$ cells, the sampling variance is estimated from within-group expression variation:

$$
\widehat\delta_{ic}=s_{ic}^{2}/n_{ic},
$$

where $s_{ic}^{2}$ is the sample variance across cells on the chosen expression scale. It enters as an estimated covariance offset. A donor represented by fewer cells can have a noisier expression mean. [CIGMA model and sampling correction](https://www.nature.com/articles/s41586-026-10577-6)

### Genetic effects: a Gaussian column and correlated cell types

Define $\beta_{jc}=\alpha_j+\eta_{jc}$. For one cell type,

$$
\alpha\sim N_L(0,\sigma_g^2 I_L/L),\qquad
\eta_c\sim N_L(0,v_c I_L/L),
$$

$$
\boxed{\beta_c\sim N_L\left(0,\frac{\sigma_g^2+v_c}{L}I_L\right).}
$$

Within this component, every SNP has the same **prior variance** in type $c$. Its realized effect can differ in size and sign. Across types, the columns share $\alpha$:

$$
\operatorname{Cov}(\beta_c,\beta_d)=\frac{\sigma_g^2}{L}I_L
\quad(c\ne d).
$$

Equivalently, one SNP's row vector satisfies

$$
b_j\sim N_C(0,\Omega_g/L),\qquad
\Omega_g=\sigma_g^2\mathbf1\mathbf1^\top+
\operatorname{diag}(v_1,\ldots,v_C).
$$

These are column and row descriptions of the **same prior**. The variance $v_c$ can differ by type; there is no SNP-specific selection indicator in this model. The general covariance version of CIGMA replaces this structured $\Omega_g$ with a full covariance matrix.

### Integrate the effects: the covariance that is fitted

For a complete donor-by-type expression matrix, let $X$ be the $N\times L$ genotype matrix and $K=XX^\top/L$ its donor-relatedness matrix. Define

$$
\Omega_e=\sigma_e^2\mathbf1\mathbf1^\top+
\operatorname{diag}(w_1,\ldots,w_C),\qquad
D=\operatorname{diag}\{\delta_{ic}\}.
$$

With expression stacked **by cell type**, the covariance is

$$
\boxed{
\Sigma_{\rm CIGMA}
=\underbrace{\Omega_g\otimes K}_{\text{genetic}}
+\underbrace{\Omega_e\otimes I_N}_{\text{residual donor}}
+\underbrace{D}_{\text{cell sampling}}.}
$$

$D$ uses the same stacking order. Genetic covariance follows donor relatedness $K$; residual donor covariance joins measurements of the same donor; cell sampling contributes a separate diagonal term. The [implementation](https://github.com/Minhui-Chen/CIGMA/blob/5813e4ae84d7b3733dfcd938fe42d12c6b30a8aa/src/cigma/fit.py) stacks by donor, giving the equivalent order $K\otimes\Omega_g+I_N\otimes\Omega_e+D$.

For multiple genotype components, such as cis and trans, the genetic term becomes $\sum_h\Omega_{g,h}\otimes K_h$, with $K_h=X_hX_h^\top/L_h$. Each component has its own SNP set and variance parameters.

### Estimation and scientific output

Haseman–Elston (HE) regression fits second moments. Let $P$ remove fixed effects, $y^*=Py$ be projected stacked expression, and $A_q=P\Sigma_qP$ be a projected covariance basis matrix from the decomposition above. Each $\Sigma_q$ multiplies one unknown variance parameter $\theta_q$. The fitting principle is

$$
\widehat\theta=\arg\min_\theta
\left\|y^*y^{*\top}-P\widehat D P-\sum_q\theta_qA_q\right\|_F^2.
$$

The unknown coefficients include $\sigma_g^2,v_c,\sigma_e^2,w_c$. Genetic and residual donor variances are estimated together, after the cell-sampling correction. The source implements this moment regression and jackknife uncertainty; REML is also available. [Fitting code](https://github.com/Minhui-Chen/CIGMA/blob/5813e4ae84d7b3733dfcd938fe42d12c6b30a8aa/src/cigma/fit.py)

For $\bar v=C^{-1}\sum_c v_c$, the genetic summaries are

$$
\text{shared variance}=\sigma_g^2,\qquad
\text{specific variance}=v_c,\qquad
\text{specificity}=\frac{\bar v}{\sigma_g^2+\bar v}.
$$

Testing $H_0:v_1=\cdots=v_C=0$ asks whether genetic effects vary across cell types for the gene. A significant result supports aggregate heterogeneity across the fitted SNP set. It does not localize that heterogeneity to a particular SNP.

## 3. CASE: infer SNP effects and their cell-type activity

### Regression and residual dependence

For $M$ cis SNPs and $N_c$ donors observed in type $c$, CASE starts from

$$
\boxed{y_c=X_c\beta_c+\varepsilon_c},\qquad
X_c\in\mathbb R^{N_c\times M},\quad
B=(\beta_1,\ldots,\beta_C)\in\mathbb R^{M\times C}.
$$

Expression and genotypes are standardized after the specified normalization and covariate adjustment. Gaussian regression residuals are independent between donors but can correlate across types measured in the same donor:

$$
\operatorname{Cov}(\varepsilon_{ic},\varepsilon_{\ell d})
=\mathbf1(i=\ell)(V_\varepsilon)_{cd}.
$$

$V_\varepsilon$ is a $C\times C$ residual covariance. CASE carries this dependence into a summary-statistic analysis through a sample-overlap adjustment. The paper's observation model and its approximation for summary statistics are distinct steps. [CASE model](https://www.nature.com/articles/s41467-026-72176-3)

### Prior on each SNP row: one pattern, several cell-type effects

Introduce a latent pattern $z_j$ for SNP $j$:

$$
z_j\sim\operatorname{Categorical}(\pi_1,\ldots,\pi_T),\qquad
\boxed{b_j\mid z_j=t\sim N_C(0,U_t)},\qquad U_T=0.
$$

SNP rows are independent conditional on the mixture parameters. The mixture $\{\pi_t,U_t\}$ is common to the SNPs in the analyzed gene. **Different SNPs can select different patterns $z_j$.** Each $U_t$ specifies the variances and correlations across cell types within one pattern.

For three cell types, an illustrative pattern is

$$
U_t=
\begin{pmatrix}
0.04&0&0.03\\
0&0&0\\
0.03&0&0.09
\end{pmatrix}.
$$

Under this pattern, SNP $j$ is active in types 1 and 3, its type-2 effect is exactly zero, and its two active effects have correlation $0.03/\sqrt{0.04\cdot0.09}=0.5$. Their sizes can differ. The all-zero pattern makes the SNP inactive in every type.

The distinction between conditional and marginal distributions is

$$
\beta_{jc}\mid z_j=t\sim N(0,(U_t)_{cc}),\qquad
\beta_{jc}\sim\sum_{t=1}^{T}\pi_tN(0,(U_t)_{cc}).
$$

Here $N(0,0)$ means a point mass at zero. Consequently,

$$
P(\beta_{jc}=0)=\sum_{t:(U_t)_{cc}=0}\pi_t.
$$

**One pattern is selected for the whole SNP row.** Its coordinates can have different variances or be zero; they do not independently select unrelated patterns. Before observing data, SNPs share the mixture law. Their posterior activity probabilities can differ once the association evidence and LD are used. [CASE prior and posterior implementation](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE_models.R)

### Summary-statistic likelihood: separate joint effects from LD

For standardized data, define the marginal association column

$$
\widehat\beta_c=X_c^\top y_c/N_c,\qquad
\widehat B=(\widehat\beta_1,\ldots,\widehat\beta_C).
$$

The LD matrix $R$ has dimension $M\times M$. It describes correlation **between SNPs**, whereas CIGMA's $K$ describes similarity **between donors**. Approximately, $E[\widehat B\mid B]=RB$.

<figure class="study-figure" style="max-width:600px">
<a href="{{base}}/assets/case-study/ld-simple.svg"><img src="{{base}}/assets/case-study/ld-simple.svg" alt="Two SNPs each affect a different cell type, but LD creates marginal associations in both types"></a>
<figcaption>LD mixes SNP effects in marginal association estimates. The example uses two SNPs and two cell types.</figcaption>
</figure>

CASE uses the asymptotic working likelihood

$$
\boxed{\widehat B\mid B\ \dot\sim\ \operatorname{MN}(RB,R,V)},\qquad
\operatorname{Cov}(\widehat B_{jc},\widehat B_{kd}\mid B)
\approx R_{jk}V_{cd}.
$$

The matrices describe different levels of variation:

| Matrix | Meaning |
|---|---|
| $U_t$ | Prior covariance of a SNP's effects, conditional on sharing pattern $t$ |
| $V_\varepsilon$ | Residual covariance in the donor-level regression |
| $V_y$ | Total phenotypic covariance across types on the standardized scale |
| $V$ | Across-type sampling covariance used in the summary likelihood |

The paper uses

$$
V_{cd}=(V_y)_{cd}\frac{N_{cd}}{N_cN_d},\qquad V_{cc}=1/N_c,
$$

where $N_{cd}$ counts donors observed in both types. This is the paper's summary-likelihood approximation; $V$ and $V_\varepsilon$ have different definitions. The approximation also assumes suitable LD for the donor groups.

### Fit the mixture, then calculate posterior support

The analysis proceeds through

$$
(\widehat B,R,\{N_c,N_{cd}\})
\longrightarrow\widehat V
\longrightarrow\{\widehat\pi_t,\widehat U_t\}
\longrightarrow P(B,z\mid\text{data})
\longrightarrow\operatorname{PIP}_{jc}.
$$

The paper estimates $V$ from weak association signals and fits mixture parameters using Monte Carlo EM. Its E-step uses MCMC for the latent SNP effects and pattern allocations. Posterior sampling with the fitted parameters gives

$$
\operatorname{PIP}_{jc}=P(\beta_{jc}\ne0\mid\widehat B,R,\widehat V,
\{\widehat\pi_t,\widehat U_t\}).
$$

LD couples SNPs in the likelihood; the row prior couples cell types. The [software interface](https://github.com/leaffur/CASE/blob/13f4fc8432ba39853128a4afc5dd8c6a151589a0/R/CASE.R) also accepts z-scores or marginal effects with standard errors. Its covariance and overlap defaults should be checked when reproducing the paper's analysis.

### From PIPs to cell-type-wise credible sets

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

## 4. Compare the priors and the questions they answer

### What is shared across SNPs, and what can vary?

| View of the effects | CIGMA shared-plus-specific model | CASE |
|---|---|---|
| Fixed type $c$, vary SNP $j$ | $\beta_{jc}\sim N(0,(\sigma_g^2+v_c)/L)$ for every SNP in the component | $\beta_{jc}$ has the same marginal mixture law across SNPs, but each SNP has its own latent pattern $z_j$ |
| Fixed SNP $j$, vary type $c$ | One Gaussian row prior with shared covariance and type-specific diagonal variances | One mixture component for the row; its covariance can select a subset of active types |
| What differs between SNPs? | Realized effect values | Realized effect values, latent sharing patterns, and posterior support |
| What is estimated? | Genetic and residual variance components | Mixture parameters and posterior SNP–type effects |

For CASE, conditioning on all SNP patterns makes the cell-type column Gaussian:

$$
\beta_c\mid z_1,\ldots,z_M
\sim N_M\!\left(0,
\operatorname{diag}\big((U_{z_1})_{cc},\ldots,(U_{z_M})_{cc}\big)\right).
$$

Some entries have zero variance, and active entries can have different variances because their SNPs selected different patterns. Averaging over the unknown patterns gives a mixture distribution for the column. This differs from CIGMA's single Gaussian column with one variance shared across its SNP entries.

A Gaussian prior does not require equal realized SNP effects. A mixture prior does not assign a separately estimated prior to every SNP. Both distinctions concern the distribution before observing the data; LD can induce posterior dependence between SNPs even when their prior effects are independent.

### A connection at the genetic-covariance level

For an aligned SNP set and genotype scale, suppose the SNP rows are independent with zero mean and common covariance $S$. The genetic component alone, $A=XB$, then satisfies

$$
\operatorname{Cov}(A_{ic},A_{\ell d}\mid X)
=\sum_j x_{ij}x_{\ell j}S_{cd}.
$$

CIGMA uses $S=\Omega_g/L$. CASE's mixture implies $S=\sum_t\pi_tU_t$. This gives a way to compare their implied genetic covariance. The full observation models still contain their own residual and sampling terms.

CIGMA's HE estimator uses second moments, so its estimation principle can extend beyond a literally Gaussian collection of causal effects. CASE's posterior activity probabilities also depend on the mixture's zero components and distributional shape. Equal second moments alone cannot determine those probabilities.

### Choose the analysis from the target

| Scientific question | Target | Relevant analysis |
|---|---|---|
| How much genetic variation is shared or heterogeneous across types? | $\sigma_g^2$, $v_c$, and specificity | CIGMA variance estimates and uncertainty |
| Does the gene show heterogeneous genetic regulation? | $H_0:v_1=\cdots=v_C=0$ | CIGMA gene-level test |
| Which variant could regulate the gene in type $c$? | $P(\beta_{jc}\ne0\mid\text{data})$ | CASE PIPs and credible sets |
| Which types support activity of a variant? | Nonzero coordinates of $b_j$ | CASE joint posterior |
| Does a type contain at least one supported regulatory signal? | A reported credible set for that gene and type | CASE eGene call |

Shared activity and equal effects are different targets. A SNP with effects $(0.1,0.3,0.7)$ is active in every type and varies in magnitude. CASE can support its activity pattern; CIGMA can detect aggregate heterogeneity. A missing discovery in either method can reflect limited precision.

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

**Revision: October 1, 2026. Data and computation.** The gene-overlap counts are transcribed from the September 28 OneK1K progress report and have not been re-estimated here. The LD and credible-set diagrams are teaching examples. [Figure-generation code]({{base}}/assets/case-study/make_diagrams.py). Source inspection refers to CIGMA commit `5813e4a` and CASE commit `13f4fc8`; no new OneK1K model fit is reported.
