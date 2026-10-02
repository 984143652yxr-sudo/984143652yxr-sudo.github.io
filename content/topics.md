# Study Topics

<p class="lead">Study notes on statistical genetics and the mathematical mechanisms of language models.</p>

These sections organize the questions in my study diary. The linked entries document learning and ongoing replication work.

<section class="research-block" markdown="1">

## Cell-type-specific genetic regulation

The same genome is present across many cell types, but the relationship between genotype and expression can depend on cellular context. I am studying how variance-component models separate shared genetic regulation from cell-type-specific effects.

- How does CIGMA distinguish shared and specific genetic variance?
- How does finite cell sampling affect pseudobulk precision?
- When does a difference in statistical significance reflect power rather than a difference in biology?

**Start here:** [Genetic variation across genes, cells, and cell states]({{base}}/diary/cell-study/) — measurement, eQTLs, genetic variance, and cell-state models.

**Fine-mapping:** [CIGMA and CASE: variance components, fine-mapping, and cellular context]({{base}}/diary/case-fine-mapping/) — statistical models, credible-set construction, overlap, and biological enrichment.

**Replication:** [Replicating CIGMA on OneK1K]({{base}}/diary/cigma-onek1k/).

</section>

<section class="research-block" markdown="1">

## Reproducibility & statistical calibration

A replication is a sequence of decisions about data, models, and evidence. I want to make those decisions visible: the checks that pass, the assumptions that remain provisional, and the differences between an updated analysis and the original study.

- Keep donor, gene, and cell-type identities aligned across every input.
- Track how filtering changes the analysis population.
- Diagnose calibration before interpreting biological enrichment.

**First case study:** the [workflow notes]({{base}}/diary/cigma-onek1k/) from the OneK1K replication.

</section>

<section class="research-block" markdown="1">

## Language models: attention and representation

Mathematical derivations and small executable examples for understanding how language models represent and process a sequence.

**Position and attention:** [Rotary position embeddings: from rotations to attention]({{base}}/diary/rotary-position-embedding/) — relative-position geometry, coordinate pairing, and numerical checks.

</section>
