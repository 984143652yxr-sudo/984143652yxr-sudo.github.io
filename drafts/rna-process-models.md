# Deferred study: RNA regulatory processes

Saved from the cell study for a separate future entry. This file is outside the published content directory.

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

