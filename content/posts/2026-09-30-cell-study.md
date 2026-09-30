---
{
  "slug": "cell-study",
  "title": "How do we understand changes in single-cell expression?",
  "date": "2026-09-30",
  "date_label": "September 30, 2026",
  "status": "Study note",
  "reading_time": "9 min read",
  "summary": "From a gene's DNA to its measured RNA count: genetic effects, cellular regulation, and models for the steps in between."
}
---

## My starting question

My starting point is CIGMA and the question of how genetic effects vary across cell types. To understand that question, I first want a picture of the whole process: how a cell produces RNA, how the experiment measures it, and which part each statistical model describes.

A useful sequence is **DNA → regulation → RNA production and processing → RNA abundance → measured count**. Genotype can influence several steps. Cell type, cell state, and the experiment also shape what we observe.

## 1. From a person to a count

A **cell** contains DNA and the machinery that reads it. A **gene** is a DNA region used to produce RNA. An **RNA transcript** is a physical molecule made from that template. In single-cell RNA sequencing, each matrix entry measures RNA assigned to **one gene in one cell**.

![From a person to a gene-expression count]({{base}}/assets/cell-study/01-what-is-measured.png)

*Figure 1. Follow the process into one highlighted cell–gene entry. Gene X, Gene Y, and the counts are illustrative. The diagram combines RNA processing and export; processing occurs primarily in the nucleus.*

The experiment captures a sample of the RNA and converts it to tagged cDNA for sequencing. A cell barcode identifies the cell, sequence matching identifies the gene, and a unique molecular identifier (UMI) helps count molecules while accounting for amplification duplicates. A count of 2 means two counted molecular identifiers assigned to that cell–gene pair. [10x: from sequencing reads to counts](https://www.10xgenomics.com/blog/how-single-cell-sequencing-data-analysis-works)

The datasets join at different levels:

| Data | One row represents | Columns contain |
|---|---|---|
| Genotypes | A donor | Allele dosages at SNP positions |
| Expression | A cell | RNA counts for genes |
| Cell metadata | A cell | Donor ID, cell type, state, batch |

For donor-level analysis, cells can be grouped by donor and cell type. Summing raw counts gives a pseudobulk count profile; averaging transformed expression gives a different summary. The choice determines the response scale and its sampling uncertainty. Additional cells improve the measurement of a donor's expression; additional donors provide more genetic samples.

## 2. Genetic differences between people

A **SNP** is a position with a single-base variant. At the same position, one donor may carry A/A, another A/G, and another G/G. Counting G gives dosages 0, 1, and 2. Each donor generally carries that inherited autosomal genotype across their cell types.

![Genotype differences and a possible regulatory mechanism]({{base}}/assets/cell-study/02-genotype-effects.png)

*Figure 2. A hypothetical allele changes protein binding and transcription. The table shows expected expression on an illustrative analysis scale. The liver-cell example shows how the effect can depend on cellular context.*

For one gene in one cell type, the familiar regression question is:

$$
\text{expression} = \text{baseline} + \beta\,\text{dosage} + \text{covariate effects} + \text{residual}.
$$

Here $\beta$ describes the association per additional counted allele. Comparing this coefficient across cell types asks whether the genetic effect changes with context. Comparing cell-type means asks about baseline expression. Correlated variants and donor differences enter the analysis through variant selection, covariates, and the covariance model.

CIGMA extends the genetic question to **shared and cell-type-specific genetic variance**, using genotype information and donor-by-cell-type expression summaries. This is the main focus of my [OneK1K replication diary]({{base}}/diary/cigma-onek1k/). Its genetic variance components describe variation associated with genetic similarity across donors. [CIGMA source](https://github.com/Minhui-Chen/CIGMA)

## 3. What can change at the same genotype?

Cells with the same inherited genotype can differ in DNA accessibility, regulatory proteins, signals, and RNA stability. Even cells of the same type fluctuate in state and transcriptional activity.

![Regulation, RNA lifetime, fluctuations, and detection]({{base}}/assets/cell-study/03-fixed-genotype-mechanisms.png)

*Figure 3. Each row follows a possible route to different RNA abundance. The bottom panel holds abundance fixed and changes detection. These are teaching mechanisms; their importance depends on the gene and experiment.*

RNA abundance reflects the balance between production and removal. With constant rates at equilibrium, expected abundance equals **production rate / degradation rate**. The observed count also depends on capture and sequencing. A larger count can therefore arise through several routes, which motivates different models and experimental measurements. Regulatory mechanisms can themselves be influenced by genetics; here I am comparing cells while holding genotype fixed.

## 4. Models for the steps in the figure

The following papers give complementary views of the process. I have shortened their notation and kept the inputs, modeling idea, and interpretation together.

### Regulatory activity: SCENIC

**Question:** Which transcription-factor programs are active in different cells?

[Aibar et al., 2017, *Nature Methods*](https://www.nature.com/articles/nmeth.4463) use random-forest regression to predict each gene's expression from transcription-factor expression. They retain candidate regulator–target groups with supporting DNA-binding motifs, then score each group's activity in each cell using expression ranks.

**Inputs:** cell expression profiles and motif-reference information. **Output:** candidate regulatory networks and cell-level activity scores. This connects the regulatory-protein panel to a statistical learning problem. The inferred links guide hypotheses; direct binding and causal effects require experimental validation. Accessibility itself can be measured with ATAC-seq, providing a complementary layer of evidence.

### Bursts and capture: a telegraph model

**Question:** Do cells differ because a gene turns on more often, stays on longer, or produces RNA faster?

[Tang et al., 2023, *Bioinformatics*](https://doi.org/10.1093/bioinformatics/btad395) model a promoter switching between OFF and ON, RNA synthesis while ON, and RNA degradation. They include cell size and capture efficiency when inferring burst kinetics from count distributions.

**Inputs:** single-cell counts and assumptions or estimates for size and capture. **Output:** kinetic parameters, with identifiability depending on those assumptions. Their likelihood and moment-based approaches account for how measurement changes the distribution used to infer biological kinetics.

<details markdown="1">
<summary>The model in a few symbols</summary>

$$
\mathrm{OFF}\underset{k_{\rm off}}{\overset{k_{\rm on}}{\rightleftharpoons}}\mathrm{ON},
\qquad C\mid M,q\sim\mathrm{Binomial}(M,q).
$$

The promoter activates at rate $k_{\rm on}$ and deactivates at $k_{\rm off}$. While active, it produces RNA at a synthesis rate; each molecule has a degradation rate. $M$ is the RNA number and $C$ its sampled count, with detection probability $q$. The basic stationary telegraph model gives a Beta–Poisson count distribution; cell size and capture alter the effective synthesis scale. Absolute time scales need additional information.

</details>

### Processing and decay: scVelo

**Question:** Is a gene's RNA abundance rising or falling as cells change state?

[Bergen et al., 2020, *Nature Biotechnology*](https://doi.org/10.1038/s41587-020-0591-3) model transcription, splicing, and degradation using **unspliced and spliced RNA**. The dynamical model fits reaction rates and latent positions along a gene's trajectory by likelihood-based expectation–maximization. [Authors' model description](https://scvelo.readthedocs.io/en/latest/about.html)

**Inputs:** separate unspliced and spliced measurements across cells spanning relevant states. **Output:** inferred expression-change directions and latent progression under the kinetic model. A gene-total count matrix alone leaves out the two RNA compartments this model uses. The inferred time scale needs external calibration for interpretation in physical time.

<details markdown="1">
<summary>The two-compartment model</summary>

$$
\frac{du}{dt}=a(t)-b u,\qquad
\frac{ds}{dt}=b u-d s.
$$

Here $u$ is unspliced RNA, $s$ is spliced RNA, $a(t)$ is transcription, $b$ is the splicing rate, and $d$ is the degradation rate. The first equation tracks newly transcribed RNA entering and leaving the unspliced compartment. The second tracks mature RNA gaining molecules through splicing and losing them through degradation. These letters are simplified from the paper's notation.

</details>

### Adding experimental time: newly synthesized RNA

**Question:** Can we measure recent production separately from RNA already present?

[Ramsköld et al., 2024, *Nature Cell Biology*](https://www.nature.com/articles/s41556-024-01486-9) developed NASC-seq2 to label newly made RNA during a defined interval. A mixture model separates labeling-related base conversions from background errors. New-RNA counts, together with degradation-rate information, support likelihood-based inference of transcriptional burst parameters.

**Inputs:** a metabolic-labeling experiment with a known labeling interval. **Output:** information about recent synthesis and burst kinetics. This adds a time reference to the observation process and helps investigate the production side of Figure 3.

## 5. What I would investigate next

| My question | Useful data and starting model |
|---|---|
| Does genotype predict expression differently across cell types? | Matched donor genotypes and expression; eQTL or variance-component model |
| Which regulatory programs differ between cells? | Expression plus motif information; SCENIC, with accessibility or perturbation data for follow-up |
| What explains variation in RNA production? | Cell counts with a capture-aware burst model; labeling data adds temporal information |
| Is expression increasing or decreasing during a transition? | Unspliced/spliced measurements and a kinetic model such as scVelo |

My next reading question is **identifiability**: when two mechanisms produce similar count distributions, what extra measurement separates them? For the OneK1K work, I will first check which RNA layers and cell-state annotations are available, then choose the model around the question and those measurements.
