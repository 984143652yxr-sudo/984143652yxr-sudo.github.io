---
{
  "slug": "cell-study",
  "title": "Modeling regulation from genotype to single-cell RNA",
  "date": "2026-09-30",
  "date_label": "September 30, 2026",
  "status": "Study note",
  "reading_time": "13 min read",
  "summary": "A statistical map of RNA measurement, shared genetic effects, regulatory kinetics, and cell-state abundance—including our GeNA study."
}
---

## Starting question

How does genetic variation change RNA expression, and how does that effect depend on the cell? I want to follow **DNA → regulation → RNA production and processing → RNA abundance → measured count**, specifying a response and a model at each step.

The page follows three questions: **what is measured** (Section 1), **how expression associates with genotype** (Section 2), and **which regulatory processes or cell states could explain the variation** (Section 3).

The equations below hold **one gene fixed** unless stated otherwise. Subscripts appear when we compare donors or cell types. This keeps the biological object and the statistical unit visible.

| Symbol | Meaning |
|---|---|
| $i$, $c$ | Donor and cell type |
| $G_i$, $Z$ | One donor's allele dosage; standardized donor-by-variant matrix |
| $M$, $\lambda$ | RNA molecules present in one cell; their expected abundance |
| $C$, $q$, $L$ | Observed gene count, detection probability, total cell UMI count |
| $X$ | Transformed expression in one cell |
| $B_{ic}$, $Y_{ic}$ | Sum of raw counts; mean transformed expression in a donor–type group |
| $S$, $W$ | Cell state; measured adjustment covariates |
| $\beta_c$, $K$ | SNP effect in type $c$; donor genetic-relatedness matrix |

## 1. From a person to an expression measurement

A **gene** is a DNA region used as a template for RNA. A **SNP** is a position with a single-base variant. A **cell** contains the molecular machinery that transcribes DNA. Repeated transcription produces multiple **RNA molecules** from a gene. Each entry of a single-cell count matrix records RNA assigned to **one gene in one cell**.

<figure class="study-figure">
<a href="{{base}}/assets/cell-study/01-what-is-measured.png"><img src="{{base}}/assets/cell-study/01-what-is-measured.png" alt="Person, cells, transcription, capture and sequencing, ending with a cell-by-gene count matrix"></a>
<figcaption>Figure 1. Follow one gene from DNA to a measured count. Illustrative molecules and counts. Click any figure to enlarge.</figcaption>
</figure>

A cell barcode identifies the cell, sequence alignment identifies the gene, and a UMI helps distinguish captured molecules from amplification duplicates. Capture, conversion to cDNA, sequencing, and assignment determine which molecules reach the matrix. RNA processing occurs primarily in the nucleus before export. [10x: reads to counts](https://www.10xgenomics.com/blog/how-single-cell-sequencing-data-analysis-works)

### Observation model: abundance and detection

For one gene in one cell, an elementary sampling model is

$$
M\mid\lambda\sim\operatorname{Poisson}(\lambda),\qquad
C\mid M,q\sim\operatorname{Binomial}(M,q).
$$

Poisson thinning gives

$$
C\mid\lambda,q\sim\operatorname{Poisson}(q\lambda),\qquad E[C]=q\lambda.
$$

Here $\lambda$ describes biology and $q$ summarizes detection. For example, $(\lambda,q)=(20,0.1)$ and $(10,0.2)$ both give an expected count of 2. Count data alone identify their product under this model. Independent detection information helps separate them.

The Poisson assumption describes a snapshot count distribution. A common extension allows extra variation:

$$
C\sim\operatorname{NB}(\mu,\phi),\qquad
E[C]=\mu,\quad \operatorname{Var}(C)=\mu+\phi\mu^2.
$$

This parameterization separates the mean from overdispersion. Section 3 introduces mechanisms that can generate heterogeneous RNA abundance.

### Data arrays and donor summaries

| Array | Rows | Columns |
|---|---|---|
| Genotypes | Donors | SNP dosages, usually 0, 1, 2; imputed dosages can be fractional |
| RNA counts | Cells | Genes |
| Metadata | Cells | Donor ID, cell type, state, batch |

For $N$ donors, $P$ variants, $H$ cells, and $Q$ genes, the array sizes are $N\times P$ and $H\times Q$. Donor IDs join the two arrays. Each donor contributes multiple cell measurements of the same inherited genotype.

For a cell with $L>0$, one familiar transformation is

$$
X=\log\left(1+10^4 C/L\right).
$$

Within donor $i$ and type $c$, let the sum run over that group's $n_{ic}$ cells:

$$
B_{ic}=\sum_{\text{cells}}C,\qquad
Y_{ic}=\frac{1}{n_{ic}}\sum_{\text{cells}}X.
$$

$B$ is a pseudobulk raw count, analyzed with an appropriate library-size adjustment. $Y$ is mean transformed expression. A group containing counts 1 and 3 has $B=4$; its $Y$ also depends on the cells' library sizes. The selected summary fixes the interpretation of a regression coefficient.

More cells improve each donor's expression measurement; more donors supply additional genetic observations.

## 2. Genotype, eQTLs, and shared regulatory effects

At a fixed autosomal SNP, A/A, A/G, and G/G correspond to $G_i=0,1,2$ when counting G. A donor generally carries that inherited genotype across cell types. A variant can lie in coding sequence, an intron, or a nearby regulatory region. An eQTL is an association between a genomic locus and expression; cis/trans terminology concerns its relationship to the target gene. Coding status and cis/trans classification describe different properties.

<figure class="study-figure">
<a href="{{base}}/assets/cell-study/02-genotype-effects.png"><img src="{{base}}/assets/cell-study/02-genotype-effects.png" alt="Same SNP position across three donors; allele-dependent transcription differs across cell types"></a>
<figcaption>Figure 2. A hypothetical SNP effect depends on cellular context. The numbers represent expected expression on an illustrative analysis scale.</figcaption>
</figure>

### One SNP: an effect on an expression response

For one gene and cell type, regress donor summaries on dosage:

$$
Y_{ic}=\mu_c+\beta_cG_i+W_i^\top\gamma_c+\varepsilon_{ic}.
$$

$\beta_c$ is the expression difference per additional counted allele, conditional on the covariates. These can include ancestry and technical factors. Related donors and repeated cell-type measurements require corresponding covariance terms.

The figure's example has

$$
E[Y_{i,\mathrm T}\mid G_i]=2+2G_i,\qquad
E[Y_{i,\mathrm L}\mid G_i]=2.
$$

The T-cell SNP effect is 2 and the liver-cell effect is 0. We can separately test average expression differences, $\beta_{\mathrm T}=0$, or $\beta_{\mathrm T}=\beta_{\mathrm L}$. An associated SNP may tag a correlated causal variant through linkage disequilibrium; functional evidence helps locate the regulatory mechanism.

### Multiple SNPs: a modeled genetic contribution

For a chosen set of variants, the genetic component is a weighted sum:

$$
a_{ic}=\sum_{j=1}^{P}Z_{ij}\beta_{jc},\qquad
\mathbf y_c=\mu_c\mathbf1+W\gamma_c+Z\boldsymbol\beta_c+\boldsymbol\varepsilon_c.
$$

The gene is the DNA feature; $a_{ic}$ is a contribution to its expression measurement. For standardized predictors $(1,-1)$ and coefficients $(0.4,0.1)$, this contribution is $0.3$ on the chosen expression scale. Predictor variants may be selected from a defined cis window.

### Shared and cell-type-specific variance

A compact random-effect construction writes each variant's effect as

$$
\beta_{jc}=\alpha_j+\eta_{jc},\qquad K=ZZ^\top/P.
$$

Assume zero-mean coefficients with $\operatorname{Var}(\alpha_j)=\sigma^2_{\rm shared}/P$ and $\operatorname{Var}(\eta_{jc})=v_c/P$. Take these coefficients independent across variants, between components, and across types for the specific component. Conditional on $Z$, this gives

$$
\operatorname{Cov}(a_{ic},a_{\ell d})
=K_{i\ell}\left[\sigma^2_{\rm shared}+\mathbf1(c=d)v_c\right].
$$

Genetically similar donors have correlated genetic contributions. Within a type, the covariance includes shared and specific components; across types, the shared component remains.

The expression model combines this genetic contribution with covariates and a residual:

$$
Y_{ic}=\mu_c+W_i^\top\gamma_c+a_{ic}+r_{ic}.
$$

Here $r_{ic}$ collects remaining biological and measurement variation. Residuals from different cell types of the same donor can be correlated; a residual covariance matrix $R$ represents that dependence. The genetic component is distinguished through the relatedness matrix $K$, which connects genetically similar donors.

This is a teaching construction for **CIGMA's shared/specific genetic-variance question**. The [CIGMA implementation](https://github.com/Minhui-Chen/CIGMA) and [replication diary]({{base}}/diary/cigma-onek1k/) give the method-specific specification.

## 3. Regulatory processes and cell-state models

At the same genotype, cells can differ in chromatin accessibility, regulatory proteins, signals, and RNA lifetime. These processes can also mediate genetic effects.

Read Figure 3 from regulation to measurement: **regulators change RNA production; processing and degradation change the RNA present; detection produces the observed count**. The first five models below describe parts of this process. The final two ask how expression and cell-state abundance associate with genotype.

<figure class="study-figure">
<a href="{{base}}/assets/cell-study/03-fixed-genotype-mechanisms.png"><img src="{{base}}/assets/cell-study/03-fixed-genotype-mechanisms.png" alt="DNA accessibility, regulatory proteins, signals and RNA degradation alter abundance; capture alters measured counts"></a>
<figcaption>Figure 3. Biological regulation changes abundance; detection changes the sampled count. Each mechanism suggests a different statistical question.</figcaption>
</figure>

### 3.1 Production and degradation: which rate changed?

**Process:** RNA is produced and removed. **Quantity modeled:** the expected number of molecules of one gene in a cell.

Let $\lambda(t)$ be expected abundance, $a$ the production rate (molecules per unit time), and $d$ the per-molecule degradation rate (per unit time):

$$
\frac{d\lambda}{dt}=a-d\lambda,\qquad
\lambda_{\rm equilibrium}=a/d,\qquad E[C]=qa/d.
$$

At equilibrium, production balances removal. Doubling $a$ doubles expected abundance; doubling $d$ halves it. Increasing detection $q$ raises the expected measured count while leaving the RNA present in the cell unchanged. A snapshot mean identifies a combination of rates. Time-resolved measurements help distinguish them.

### 3.2 Transcriptional bursts: a switching promoter

**Process:** a promoter, the DNA region controlling transcription initiation, switches between inactive and active states. **Data:** a distribution of counts across cells.

The telegraph model represents that switching and the resulting molecule counts:

$$
\mathrm{OFF}\underset{k_{\rm off}}{\overset{k_{\rm on}}{\rightleftharpoons}}\mathrm{ON},
\qquad
M\longrightarrow M+1\ \text{ at rate }a\text{ while ON},
\qquad
M\longrightarrow M-1\ \text{ at rate }dM.
$$

At equilibrium,

$$
E[M]=\frac{a}{d}\frac{k_{\rm on}}{k_{\rm on}+k_{\rm off}},\qquad
C\mid M,q\sim\operatorname{Binomial}(M,q).
$$

$k_{\rm on}$ and $k_{\rm off}$ are switching rates. The ratio $k_{\rm on}/(k_{\rm on}+k_{\rm off})$ is the long-run fraction of time active; $a$ controls production while active; $d$ controls RNA lifetime. Counts across cells inform these parameters through their distribution. [Tang et al. (2023)](https://doi.org/10.1093/bioinformatics/btad395) incorporate cell size and capture into burst-kinetic inference using likelihood and moment methods. Their model motivates adjusting the effective synthesis scale for both. Absolute rates require a time reference or externally specified degradation scale.

### 3.3 Regulatory proteins: SCENIC

**Process:** transcription factors bind regulatory DNA and influence target-gene transcription. **Data:** expression across cells plus DNA-motif reference information.

For one target gene, the prediction step can be represented schematically as

$$
X=f(T_1,\ldots,T_p)+\epsilon,
$$

where $T_1,\ldots,T_p$ are transcription-factor expression measurements. [SCENIC, Aibar et al. (2017)](https://www.nature.com/articles/nmeth.4463) uses random-forest prediction, motif enrichment to support candidate regulator–target groups, and AUCell expression-rank scores for each cell's regulatory program activity.

The statistical target is a candidate regulatory network and its activity across cells. Motif evidence supports sequence compatibility; ATAC-seq adds evidence about accessible DNA. Binding and perturbation experiments help test the proposed regulation. The displayed equation summarizes the prediction task.

### 3.4 RNA processing and changing state: scVelo

**Process:** newly transcribed RNA is spliced into mature RNA, which later degrades. **Data:** separate unspliced and spliced measurements from cells spanning a transition.

Let $u$ and $s$ denote unspliced and spliced RNA abundance. The two-compartment model is

$$
\frac{du}{dt}=a(t)-bu,\qquad
\frac{ds}{dt}=bu-ds.
$$

Here $b$ is the splicing rate and $d$ the degradation rate. Transcription adds unspliced RNA at rate $a(t)$; splicing transfers it at rate $bu$; mature RNA decays at rate $ds$. Thus $bu-ds$ describes the modeled direction of mature-RNA change.

[scVelo, Bergen et al. (2020)](https://doi.org/10.1038/s41587-020-0591-3) fits reaction rates and latent cell times/states through likelihood-based expectation–maximization. It uses cells spanning a trajectory and distinct spliced/unspliced layers. The inferred direction depends on the kinetic model, and physical time needs calibration. [Authors' equations and fitting procedure](https://scvelo.readthedocs.io/en/latest/about.html)

### 3.5 A measured time interval: new RNA

**Process:** label RNA produced during a known interval and measure how much new RNA accumulates.

[NASC-seq2, Ramsköld et al. (2024)](https://www.nature.com/articles/s41556-024-01486-9) labels newly synthesized RNA during a known interval. A mixture model distinguishes labeling-related base conversions from background errors; new-RNA counts and degradation information support burst inference.

To see the information gained, a constant-rate teaching model starting with zero labeled molecules gives

$$
E[M_{\rm new}(t)]=\frac{a}{d}(1-e^{-dt}).
$$

This follows from the production–decay equation. The known interval $t$ supplies a time reference. This equation explains the experimental idea; the paper's inference additionally models labeling and transcriptional bursting.

### 3.6 Expression within a cell state: a genotype-by-state model

**Question:** does a SNP have a stronger expression effect in activated cells? **Data:** gene counts, donor genotypes, and cell-state measurements.

A cell type is a broad identity, such as a T cell. A continuous state $S$ can describe activation or progression within that type. For one gene, a teaching count model for cell $k$ from donor $i$ is

$$
C_{ik}\sim\operatorname{NB}(\mu_{ik},\phi),
$$

$$
\log\mu_{ik}=\log L_{ik}+\alpha+\beta G_i+f(S_{ik})
+G_i h(S_{ik})+W_{ik}^{\top}\gamma+u_i.
$$

Read the mean model term by term:

| Term | Role |
|---|---|
| $\log L$ | Adjust for the cell's library size, with $L>0$ |
| $\beta G$ | Baseline genotype association |
| $f(S)$ | Expression changes with state |
| $G h(S)$ | The genotype effect changes with state |
| $W^\top\gamma$ | Adjust for measured covariates |
| $u_i\sim N(0,\tau^2)$ | Account for cells sharing a donor |

The count distribution is conditional on these predictors and the donor effect. At state $s$, the genetic effect on the log-mean scale is $\beta+h(s)$. For $h(s)=\theta s$, moving from state 0 to state 1 changes that effect from $\beta$ to $\beta+\theta$.

The response is **one gene's expression conditional on state**. State-dependent eQTL studies motivate this question; see [Nathan et al. (2022)](https://www.nature.com/articles/s41586-022-04713-1). The equation above is an adaptable negative-binomial formulation. State estimates derived from expression require care about target-gene leakage and uncertainty; donor-level validation and sensitivity analyses help assess this.

### 3.7 GeNA: genotype and the abundance of cell states

**Question:** do donors with different genotypes have different proportions of cells in particular states?

Our **OneK1K → GeNA study notebook** models each donor's distribution of cells across transcriptional neighborhoods. A neighborhood is a local region of a cell-state graph. The sequence is **cell expression → cell graph → donor neighborhood abundances → abundance PCs → genotype association**. [Rumker et al. (2024)](https://doi.org/10.1038/s41588-024-01909-1)

<figure class="study-figure study-figure-wide">
<a href="{{base}}/assets/GeNA_understanding.webp"><img src="{{base}}/assets/GeNA_understanding.webp" alt="GeNA workflow from cell-state neighborhoods to donor abundance principal components and genetic association"></a>
<figcaption>Figure 4. GeNA overview. Source: Rumker et al., Nature Genetics (2024), Figure 1; open the linked paper for its full caption.</figcaption>
</figure>

**Step 1 — describe each donor's cell distribution.** Let $w_{km}$ be the soft contribution of cell $k$ to neighborhood $m$. A compact representation of the donor-by-neighborhood abundance matrix is

$$
Q_{im}=\frac{1}{n_i}\sum_{k\in\text{donor }i}w_{km},\qquad
\widetilde Q=UDV^\top.
$$

The entry $Q_{im}$ summarizes donor $i$'s relative abundance around neighborhood $m$. Our tutorial constructs the weights through graph diffusion.

**Step 2 — summarize abundance patterns.** Neighborhood QC and covariate/batch adjustment produce $\widetilde Q$. Its SVD gives donor PC vectors in $U$: each column represents a pattern of abundance differences across donors.

**Step 3 — test genotype association.** For one SNP, regress each retained PC on dosage:

$$
U_{ir}=a_r+b_rG_i+e_{ir},\qquad
T_r=\widehat b_r/\operatorname{se}(\widehat b_r),
$$

$$
X(k)=\sum_{r=1}^{k}T_r^2,\qquad
p(k)\approx\Pr\{\chi_k^2\ge X(k)\}.
$$

$b_r$ is the association with abundance PC $r$. Summing the squared standardized coefficients combines evidence across the first $k$ PCs. The joint test asks whether genotype associates with the abundance profile. Our inspected implementation adjusts the NAM before PCA and fits these forward regressions without extra covariates. A sensitivity analysis can explicitly adjust both genotype and phenotype. Calibration depends on the design; PC orthogonality alone supplies limited distributional guarantees.

For $J$ prespecified candidate values of $k$, the source's selection correction is

$$
p_{\rm GeNA}=1-\left(1-\min_k p(k)\right)^J.
$$

The nested tests share PCs, so this analytical procedure benefits from donor-level null checks. A neighborhood map shows where abundance increases or decreases with genotype. The denominator matters: within-NK analysis concerns relative NK states, while a PBMC-wide analysis concerns the wider mixture. [GeNA source and joint test](https://github.com/immunogenomics/GeNA)
