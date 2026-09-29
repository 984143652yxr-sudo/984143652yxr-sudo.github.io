---
{
  "slug": "cigma-onek1k",
  "title": "Replicating CIGMA on OneK1K",
  "date": "2026-09-29",
  "date_label": "September 29, 2026",
  "status": "Working draft",
  "reading_time": "7 min read",
  "summary": "From single-cell counts to shared and cell-type-specific genetic variance: a record of the workflow, debugging lessons, and questions still open."
}
---

<div class="notice"><strong>About this entry.</strong> This draft brings together the local tutorial notebooks and September 28 progress report. Numbers below are reported in that progress report; they have not been independently recomputed for this webpage. The inspected notebooks contain no saved execution outputs. This is a study record, not a claim of a fully validated replication.</div>

## The question that started this

When genetic regulation looks different across immune cell types, how much of that difference is biological, and how much reflects the precision of our measurements?

My starting point is [CIGMA](https://github.com/Minhui-Chen/CIGMA), a method for decomposing shared and cell-type-specific genetic effects on gene expression. I am using the OneK1K analysis as a way to learn the method from its inputs through to its interpretation.

> The useful output of a replication is more than a number: it is an explanation of how that number was produced.

[TOC]

## 1. Understand the units before the model

The notebook starts by separating four things that are easy to conflate: a donor, a cell, a cell type, and a gene. Cells are nested within donors. The genetic analysis needs donor-level alignment even though expression is initially observed cell by cell.

| Object | Role in the workflow | Check to preserve |
|---|---|---|
| Donor | Independent person with expression and genotype data | The same donor order in every input |
| Cell | One measured transcriptome | Correct donor and cell-type labels |
| Donor × cell type | Unit of pseudobulk expression | Sufficient cells and recorded sampling uncertainty |
| Gene | One outcome to fit | Stable gene ID and matching cis region |

The tutorial is organized into a metadata pass, a small pilot, and a production run. That sequence makes sense to me: a small pilot is where the meaning of each array should become explicit.

## 2. Freeze the cohort and cell-type mapping

The progress report describes an updated OneK1K expression release with finer cell-type labels. Those labels are mapped to seven analysis groups: `CD4_NC`, `CD4_ET`, `CD8_NC`, `CD8_ET`, `NK`, `B_IN`, and `B_MEM`.

This mapping is provisional and deserves biological review. It is one of the reasons this is an updated-data replication rather than an exact rerun of every original input.

The report records **785 retained donors** after requiring more than 10 cells in every selected type. It identifies memory B cells as an important constraint on retention. A filtering rule therefore changes both measurement precision and the population being analyzed.

**Next check:** make the donor attrition table part of the reproducibility record and compare reasonable alternatives to the cell-count threshold.

## 3. Build pseudobulk and preserve its uncertainty

The local workflow transforms each cell’s raw counts with `log10(CP10K + 1)`, then computes the mean within each donor and cell type. It also estimates the variance of that mean from the cells contributing to it.

```text
Y[i, c]     = mean of transformed expression in donor i, type c
delta[i, c] = within-group sample variance / number of cells
```

The key distinction in my notes is between uncertainty from sampling a finite number of cells and residual variation between donors. They are different parts of the model.

The tutorial also scales expression using overall pseudobulk. When expression is divided by a standard deviation, its sampling variance must be divided by the corresponding variance. Keeping these transformations paired is a useful implementation check.

## 4. Align expression, genotypes, and gene coordinates

This is the stage where an apparently successful computation can still be answering the wrong question. The local workflow combines expression gene identifiers with GRCh37 genotype coordinates, so it uses a matching gene annotation joined by stable Ensembl ID.

Each gene needs its own cis genotype matrix and derived relatedness matrix. The following is a **schematic check**, not a replacement for the notebook’s implementation:

```python
# Check labelled objects before converting to NumPy arrays.
assert expression_donors.equals(genotype_donors)
assert expression_donors.equals(covariate_donors)
assert expression_gene_id == genotype_window_gene_id
assert expression_cell_types.equals(proportion_cell_types)
```

The question I want to be able to answer at any point is: “What biological entity does this row or column represent?”

## 5. Fit the model, then inspect the evidence

The report describes a gene-wise CIGMA-HE analysis with donor jackknife inference. Its covariates include sex, age groups, overall-pseudobulk PC1, genotype PCs, and sequencing pool.

For interpretation, I keep three distinctions visible:

- **Shared versus specific genetic variance:** common regulation across the selected types versus additional type-specific components.
- **Genetic versus residual specificity:** different model components, with different interpretations.
- **Specific genetic regulation versus differential expression:** variation in genetic effects and variation in average expression ask different questions.

I also keep raw variance-component estimates. An unconstrained method-of-moments estimate can be negative; truncating every such estimate to zero would change the summary. Ratios need separately stated denominator checks.

## What the debugging taught me

The September 28 report records several concrete failures and fixes. These are worth keeping in the diary because they explain why the validation steps exist.

| Recorded problem | Lesson for the next run |
|---|---|
| Variant-ID shell escaping caused different alleles to share an ID | Inspect a few generated IDs before running duplicate removal |
| The chunk runner did not parse the PSAM `#IID` header | Validate sample-file parsing in the pilot |
| A pilot paired one gene’s expression with another gene’s genotype matrix | Tie every fit to a checked stable gene ID |
| A covariate was missing and a PC job failed on relative paths | Audit the completed covariate matrix and use explicit input paths |

These checks belong before the full run. They are not merely cosmetic cleanup after the results arrive.

## A preliminary checkpoint

The figures below are transcribed from the local progress report. The underlying result tables and job logs still need to be linked and audited before treating this as a verified result release.

| Quantity | Reported checkpoint |
|---|---:|
| Retained donors | 785 |
| Analysis cell types | 7 |
| Genes passing the expression/annotation stage | 10,570 |
| Genes fitted | 10,561 |
| Genes without cis SNPs | 9 |
| cs-eGenes at the reported Bonferroni threshold | 159 |
| cs-eGenes at the reported 5% FDR threshold | 347 |

I am keeping this checkpoint descriptive. Similarity to a published headline count does not establish calibration or equivalence of the analysis. Differences in donors, annotations, cell-type mapping, and preprocessing need to be documented first.

## What I want to check next

1. **Freeze the provenance.** Record data release identifiers, the CIGMA commit, environment versions, the final notebook, and result-table checksums.
2. **Audit the testing universe.** Account for every eligible gene and state the multiple-testing family explicitly, including how genes without usable cis SNPs are handled.
3. **Examine calibration.** Run and inspect the planned negative controls, including donor-label permutation of the genetic relatedness matrix.
4. **Test sensitivity to cell abundance.** Assess whether cohort retention and specificity calls change under defensible sampling and filtering choices.
5. **Interpret biology after the audit.** The downstream notebook outlines comparisons with expression differences, gene features, and controls matched on detection-related information. Those remain analyses to validate, not established conclusions here.

## Sources and notebook record

- **Method:** [CIGMA source repository](https://github.com/Minhui-Chen/CIGMA).
- **Paper:** [Cell-type-specific eQTLs underlie the genetic architecture of complex traits](https://doi.org/10.1038/s41586-026-10577-6).
- **Original cohort analysis code:** [OneK1K phase 1](https://github.com/powellgenomicslab/onek1k_phase1).
- **Local working files consulted:** `OneK1K_CIGMA_tutorial.ipynb`, `OneK1K_CIGMA_tutorial_3.ipynb`, `OneK1K_CIGMA_downstream.ipynb`, and `cigma_onek1k_progress.tex`. These files are not bundled with this webpage. The final notebook version has not yet been selected for publication.

<details markdown="1">
<summary>How I plan to keep this diary</summary>

For each new entry: start with one question, record what I tried, keep the useful code or figure, explain what changed my understanding, and end with the next check. Distinguish intended analyses from completed runs and preliminary observations from validated conclusions.

</details>
