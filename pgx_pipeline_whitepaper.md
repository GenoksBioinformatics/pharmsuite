# A Force-Calling, Guideline-Concordant Pharmacogenomics Pipeline on the Ilyome Platform: Design, Implementation, and Interpretive Reporting

**Technical Whitepaper**

*Ilyome Bioinformatics Platform — Pharmacogenomics (PGx) Module*

---

## Abstract

Pharmacogenomics (PGx) translates an individual's germline genotype into actionable guidance on drug selection and dosing, yet its clinical utility depends critically on the completeness and reproducibility of the underlying genotyping. A recurrent failure mode in sequencing-based PGx is the *silent no-call*: when a variant-defining position is simply absent from a variant call format (VCF) file—because no alternate allele was observed and the position was therefore never emitted—downstream star-allele (haplotype) matchers cannot distinguish a true reference call from missing data, and may default to a lower-confidence or incorrect diplotype. Here we describe the pharmacogenomics module of the Ilyome platform, which addresses this failure mode by combining hardware-accelerated alignment (Illumina DRAGEN v4.3) with position-targeted **force-calling** using GATK HaplotypeCaller (v4.6.1) across the complete set of PGx allele-definition positions on the GRCh38 reference assembly. The resulting fully-populated VCF is normalized by the PharmCAT preprocessor and interpreted by PharmCAT (v3.2.0), yielding, for each of 22 pharmacogenes, a diplotype call, allele-level functional assignment, a predicted metabolizer phenotype, and Clinical Pharmacogenetics Implementation Consortium (CPIC)–concordant prescribing guidance. The panel comprises **1,206 distinct variant-defining genomic positions across the 22 panel genes** (1,207 records in the underlying definition file, one of which carries no resolved gene annotation). We detail the pipeline architecture, the rationale for force-calling, the interpretive logic, and an analytical-validation framework based on GeT-RM/1000 Genomes reference materials. Current limitations—most notably SNV-only resolution of the structurally complex *CYP2D6* locus—are stated explicitly, together with the planned integration path.

**Keywords:** pharmacogenomics, star alleles, PharmCAT, force-calling, CPIC, GRCh38, clinical decision support, DRAGEN

---

## 1. Background

### 1.1 The clinical rationale for pharmacogenomics

Interindividual variability in drug response is substantial, and a large fraction of it is heritable. The overwhelming majority of individuals carry at least one clinically actionable pharmacogenetic variant, meaning that PGx is not a rare-disease concern but a population-scale one. The Clinical Pharmacogenetics Implementation Consortium (CPIC), founded in 2009 as a shared project of the Pharmacogenomics Research Network and the Pharmacogenomics Knowledgebase (PharmGKB), was established specifically to remove the principal barrier to clinical adoption: the difficulty of translating a laboratory genotype into an actionable prescribing decision. CPIC now curates evidence-based, peer-reviewed gene–drug guidelines spanning dozens of genes and well over one hundred drugs, and these guidelines have become the de facto international standard for clinical PGx interpretation.

### 1.2 Star alleles, diplotypes, and the genotype-to-phenotype cascade

Pharmacogene function is conventionally described using **star (\*) allele nomenclature**, standardized by the Pharmacogene Variation (PharmVar) Consortium. A star allele is a haplotype defined by one or more variant positions; the pair of star alleles carried by an individual constitutes the **diplotype** (e.g., `CYP2C19 *1/*2`). Each star allele is assigned a **clinical function** (e.g., *normal*, *decreased*, *no function*, *increased*), and the combination of the two alleles' functions yields a **predicted metabolizer phenotype** (e.g., Normal, Intermediate, or Poor Metabolizer) or, for genes handled on an activity-score basis, a numeric activity value. This phenotype is then mapped to guideline recommendations. The interpretive cascade is therefore:

> **variant genotypes → star-allele (diplotype) call → allele functional assignment → predicted phenotype → drug-level recommendation**

Correctly executing the first step is a hard prerequisite for every subsequent one.

### 1.3 The missing-position problem in sequencing-based PGx

Star-allele matchers such as PharmCAT infer haplotypes by examining the genotype at every position that participates in an allele definition. Standard variant-calling workflows are, by design, *variant*-centric: a position at which the sample is homozygous reference and no alternate allele is observed is typically not written to the VCF at all. From the matcher's perspective this creates an irreducible ambiguity—an absent record could mean "confidently reference" or "not interrogated / no data." Because a missing position can change both the diplotype and the resulting phenotype, the conservative behavior is to down-weight or withhold the call. In panels where allele definitions depend on many positions (see §4), even a small number of unemitted reference sites can degrade an otherwise correct genotype into a no-call.

The Ilyome PGx module is built around eliminating this ambiguity at the source.

---

## 2. Pipeline overview

The Ilyome PGx module is a linear, auditable pipeline composed of four stages:

1. **Alignment and primary variant context (DRAGEN v4.3):** reads are mapped to GRCh38 and a per-sample CRAM is produced.
2. **Targeted force-calling (GATK HaplotypeCaller v4.6.1):** every PGx allele-definition position is genotyped explicitly, so that both variant and reference-homozygous calls are emitted.
3. **VCF preprocessing (PharmCAT preprocessor):** the force-called VCF is normalized, split, and aligned to the PharmCAT position specification.
4. **Interpretation and reporting (PharmCAT v3.2.0):** diplotype calling, functional/phenotype assignment, and guideline-concordant report generation.

A schematic of the data flow:

```
  Sequencing reads
        │
        ▼
┌──────────────────────┐
│  DRAGEN v4.3          │   Mapping/alignment to GRCh38
│  (alignment → CRAM)   │
└──────────┬───────────┘
           │  CRAM (GRCh38)
           ▼
┌──────────────────────────────────┐
│  GATK HaplotypeCaller v4.6.1      │   Force-call across the union of
│  (position-targeted force-call)   │   PGx allele-definition positions
└──────────┬───────────────────────┘   → variant AND reference genotypes emitted
           │  fully-populated VCF (GRCh38)
           ▼
┌──────────────────────┐
│  PharmCAT preprocessor│   Normalize, left-align, split multiallelics,
│                       │   restrict to PharmCAT positions.vcf
└──────────┬───────────┘
           │  PharmCAT-ready VCF
           ▼
┌────────────────────────────────────────────┐
│  PharmCAT v3.2.0                            │
│   ├─ Named Allele Matcher → diplotype       │
│   ├─ Phenotyper → allele function + phenotype│
│   └─ Reporter → CPIC/DPWG recommendations    │
└──────────┬─────────────────────────────────┘
           │
           ▼
   Ilyome PGx clinical report
```

---

## 3. Materials and methods

### 3.1 Reference assembly and alignment

All analyses are performed against the **GRCh38** human reference assembly, which is the assembly expected by PharmCAT's position and allele-definition files. Read alignment and CRAM generation are performed with **Illumina DRAGEN v4.3**, a hardware-accelerated secondary-analysis platform whose germline algorithms—built on multigenome (pangenome) mapping and machine-learning-based filtering—have been independently benchmarked at scale, achieving high per-class accuracy across single-nucleotide variants, indels, short tandem repeats, structural variants, and copy-number variants. Using DRAGEN for the alignment/context layer provides a fast, reproducible, and accurately-mapped substrate for the targeted genotyping stage that follows.

### 3.2 Targeted force-calling with GATK HaplotypeCaller

The core design decision of the pipeline is to **force a genotype at every PGx allele-definition position**, rather than relying on a variant-only VCF. This is implemented with **GATK HaplotypeCaller v4.6.1**, GATK's assembly-based caller, which reconstructs candidate haplotypes via a de Bruijn-like graph over each active region and assigns genotype likelihoods with a pair-HMM read-likelihood model—an approach with well-characterized indel accuracy that is important given the presence of insertion/deletion-defining alleles in several pharmacogenes.

Force-calling is restricted to, and driven by, the PGx position list: the caller is supplied with the PharmCAT allele-definition positions as the set of alleles to genotype and as the interval list, so that **each of these positions receives an explicit call—alternate genotype where an alternate allele is supported, and reference-homozygous (0/0) genotype otherwise.** The critical consequence is that homozygous-reference PGx sites are *represented in the VCF as reference calls* rather than being silently omitted; the Named Allele Matcher can then score them as confident reference observations rather than missing data.

A representative invocation (parameters shown for transparency; the production configuration is version-pinned in the Ilyome pipeline manifest):

```bash
gatk HaplotypeCaller \
    -R GRCh38.fa \
    -I sample.cram \
    -L pharmcat_positions.vcf \
    --alleles pharmcat_positions.vcf \
    --force-call-filtered-alleles true \
    -O sample.pgx.forcecalled.vcf.gz
```

> *Note to reviewers/authors:* the exact flag set that guarantees emission of reference genotypes at every requested position (e.g., use of `--alleles`/force-call mode versus a reference-confidence/`-ERC` strategy followed by joint genotyping with non-variant sites retained) should be reproduced verbatim from the production manifest before publication, to ensure bit-for-bit reproducibility.

### 3.3 VCF preprocessing (PharmCAT preprocessor)

The force-called VCF is passed through the **PharmCAT VCF preprocessor**, the tool's recommended input-preparation stage. The preprocessor normalizes representation (left-alignment and parsimony), splits multiallelic records, reconciles the VCF against PharmCAT's `positions.vcf` specification, and flags any positions that remain missing. Running the preprocessor on a force-called VCF is synergistic: force-calling ensures the positions are present, and the preprocessor guarantees they are represented in exactly the form the Named Allele Matcher expects. This division of labor—explicit genotyping upstream, canonical normalization at the interface—is what makes the diplotype calls both complete and standards-conformant.

### 3.4 Diplotype calling, phenotype assignment, and reporting (PharmCAT v3.2.0)

Interpretation is performed by **PharmCAT v3.2.0**, an open-source pharmacogenomics clinical annotation tool maintained under the ClinPGx umbrella (Stanford University and the University of Pennsylvania). PharmCAT operates in three functional stages:

- **Named Allele Matcher** — infers, per gene, the diplotype most consistent with the observed genotypes across that gene's allele-definition positions (e.g., `CYP2B6 *1/*1`).
- **Phenotyper** — assigns each called allele its CPIC-defined clinical function and combines the two allele functions into a predicted metabolizer phenotype or activity score. For activity-score genes (e.g., *DPYD*), the output is a per-allele functional value and a summed activity value rather than a metabolizer category.
- **Reporter** — links the resulting phenotype to guideline-level annotations and prescribing recommendations (CPIC, with DPWG annotations where available), producing the human-readable report.

Two representative rows from the Phenotyper/Reporter output illustrate the two reporting modes:

| Gene | Diplotype | Allele functions | Phenotype / activity |
|------|-----------|------------------|----------------------|
| *CYP2B6* | \*1/\*1 | Two Normal-function alleles | Normal Metabolizer |
| *DPYD* | Reference/Reference | Reference (Normal function) | Activity value 1.0 → see drug section |

Here *CYP2B6* is reported as a categorical metabolizer phenotype, whereas *DPYD* is reported on the activity-value basis CPIC uses for fluoropyrimidine dosing—an important distinction, because it determines how the downstream recommendation is expressed.

### 3.5 Panel content and allele-definition coverage

The module interrogates the **22-gene** PGx panel enumerated in Table 1. The underlying definition file contains **1,207 records**, of which **1,206 map to the 22 panel genes**; the remaining record (`rs12777823`, chr10:94,645,745, G>A) carries no resolved gene annotation in the source file and is excluded from the per-gene counts below.

Each definition record corresponds to a **single, distinct genomic position** (chromosome + coordinate + reference/alternate allele): within this panel there are no positions that are represented by more than one record, so the number of records per gene *is* the number of distinct force-called positions for that gene. (An earlier internal tally of "unique rsIDs" per gene undercounted several genes — e.g., G6PD, RYR1, CYP2D6 — because a large fraction of positions in this panel, particularly indels and rarer variants, have no dbSNP rsID and are recorded with a placeholder `.`; deduplicating on rsID collapses all of a gene's rsID-less positions into one, which is a database artifact and does not reflect the number of positions actually genotyped. Table 1 below reports the correct, position-based count.) All 1,206 panel positions (plus the one unassigned record) are force-called and preprocessed, guaranteeing complete positional coverage for the Named Allele Matcher.

**Table 1.** Genes in the Ilyome PGx panel, principal CPIC-guided drug associations, and force-call target size (distinct variant-defining genomic positions per gene).

| Gene | Principal drug/therapeutic association(s) | Force-called positions |
|------|-----------------------------------------|:----------------------:|
| *RYR1* | Volatile anesthetics, succinylcholine (malignant hyperthermia) | 313 |
| *G6PD* | Rasburicase and oxidative-stress drugs | 171 |
| *CYP2D6* | Codeine/opioids, tamoxifen, atomoxetine, antidepressants † | 156 |
| *CFTR* | Ivacaftor | 93 |
| *DPYD* | Fluoropyrimidines (5-FU, capecitabine) | 83 |
| *CYP2C9* | Warfarin, NSAIDs, phenytoin | 83 |
| *CYP2B6* | Efavirenz, methadone | 48 |
| *NAT2* | Isoniazid, hydralazine | 47 |
| *TPMT* | Thiopurines (azathioprine, mercaptopurine) | 45 |
| *CYP3A4* | Substrate-dependent (emerging) | 43 |
| *SLCO1B1* | Statins (simvastatin) | 35 |
| *CYP2C19* | Clopidogrel, voriconazole, SSRIs, PPIs | 35 |
| *CYP4F2* | Warfarin | 20 |
| *NUDT15* | Thiopurines | 18 |
| *CYP3A5* | Tacrolimus | 5 |
| *UGT1A1* | Atazanavir, irinotecan | 4 |
| *CACNA1S* | Volatile anesthetics, succinylcholine (malignant hyperthermia) | 2 |
| *VKORC1* | Warfarin | 1 |
| *IFNL3* | Peginterferon (HCV) | 1 |
| *F5* | Factor V Leiden (thrombophilia annotation) | 1 |
| *F2* | Prothrombin/Factor II (thrombophilia annotation) | 1 |
| *ABCG2* | Rosuvastatin, allopurinol | 1 |
| **Total (22 panel genes)** | | **1,206** |
| Unassigned‡ | — | 1 |
| **Total (source file)** | | **1,207** |

† *CYP2D6* is currently resolved at SNV resolution only; see §6 (Limitations).
‡ One record (`rs12777823`, chr10:94,645,745, G>A) is present in the source definition file without a resolved gene label and is force-called but not attributed to any of the 22 panel genes above.

---

## 4. Interpretive output

For each sample, the module emits a structured report in which every panel gene is represented by a four-part interpretive record: (i) the **diplotype**, (ii) the **allele-level functional assignment**, (iii) the **predicted phenotype or activity value**, and (iv) the associated **guideline-level recommendation** (or a pointer to the relevant drug section). The report distinguishes categorical-phenotype genes from activity-value genes, preserves the provenance of each call, and—by virtue of the force-calling design—explicitly reports reference (wild-type) diplotypes such as `*1/*1` or `Reference/Reference` as confident positive findings rather than as absent data. This is significant clinically, because a confident `*1/*1` Normal Metabolizer result carries different weight than an inconclusive/no-call result for the same gene.

Because the interpretive content is generated by PharmCAT against the current CPIC guideline set, the recommendations remain traceable to their primary sources and are updated as the underlying guideline data and PharmCAT allele definitions are versioned.

---

## 5. Analytical validation framework

A PGx pipeline intended for clinical use requires demonstrated concordance against orthogonally characterized reference materials. The recommended validation design for the Ilyome PGx module, consistent with how PharmCAT itself was validated, is:

- **Reference materials.** The CDC-coordinated **Genetic Testing Reference Materials Coordination Program (GeT-RM)** Coriell samples, for which consensus pharmacogene genotypes have been established across multiple laboratories, together with corresponding **1000 Genomes Project** sequence data. GeT-RM characterization exists for the core metabolizer genes and, importantly, for *CYP2D6*.
- **Metrics.** Concordance is assessed at two levels: **diplotype-level concordance** (called diplotype vs. reference diplotype) and **phenotype-level concordance** (predicted metabolizer phenotype/activity value vs. reference). Discordances are further partitioned into those attributable to missing positions (expected to be near-zero under force-calling), to allele-definition version differences, and to true genotyping error.
- **Acceptance criteria.** A common threshold in the PGx literature is ≥95% sensitivity and specificity at the diplotype/phenotype level for genes in scope, evaluated per gene.

> *Placeholder:* concrete concordance figures from the in-house Ilyome validation run (per-gene diplotype and phenotype concordance against GeT-RM/1000 Genomes) should be inserted here as **Table 2** prior to external release. The force-calling design is expected to specifically improve concordance on genes whose definitions include many positions that are frequently homozygous-reference (and therefore prone to omission in variant-only VCFs).

---

## 6. Limitations

The pipeline force-calls SNVs and indels at the PGx allele-definition positions and interprets them with PharmCAT. Consequently, its principal limitations concern **allele classes that are not defined by, or not recoverable from, SNV/indel genotypes in a VCF**—chiefly structural variation (SV) and copy-number variation (CNV). Because PharmCAT's Named Allele Matcher does not model gene copy number and ignores structural variants in the VCF (issuing a warning when it encounters them), these classes must be resolved by dedicated callers. The subsections below state exactly which genes and which alleles are affected, so that these calls can be treated as *provisional* and cross-checked (§6.3) before a phenotype is finalized.

**6.1 *CYP2D6*: structural/CNV alleles are not resolved (SNV-only, PharmCAT research mode).** *CYP2D6*—responsible for the metabolism of roughly one-fifth of clinically used drugs—is the single most consequential limitation. It is genotyped here using PharmCAT's research-mode calling, which infers diplotypes from **SNPs/indels only**; PharmCAT itself explicitly discourages calling *CYP2D6* from a VCF because SV and CNV are central to correct phenotype prediction and cannot be represented in a VCF. Specifically, the following alleles/configurations are **not** resolved by SNV-only calling and require review:

- **Whole-gene deletion — `*5`** (no-function). The dominant hazard is *hemizygosity misrepresentation*: when `*5` is present on one chromosome and a variant allele on the other, the surviving allele's variants appear falsely homozygous in the VCF. Per PharmCAT's own worked example, a true `*5/*29` is mis-detected as `*29/*29`, changing the diplotype and potentially the phenotype.
- **Gene duplications/multiplications — `xN`** (e.g., `*1x2`, `*2x2`, `*4x2`, `*10x2`, `*17x2`, `*36x2`). Duplications have been reported on normal-function (`*1`, `*2`), decreased-function (`*10`, `*17`) and no-function (`*4`, `*36`) backbones; only `*1`, `*2`, `*4`, and `*41` have been described with ≥3 copies. Crucially, the **Ultrarapid Metabolizer** phenotype (activity score > 2.25) is driven by duplication of functional alleles and is therefore *undetectable* without CNV information—an UM sample may be silently reported as a Normal Metabolizer.
- **Tandem / *CYP2D6–CYP2D7* hybrid alleles** (e.g., `*13`, `*36`, `*68`, and tandem arrangements such as `*68+*4` or `*36+*10`), which arise from recombination between the gene and its pseudogene and cannot be reconstructed from SNV zygosity alone.

**6.2 Other loci with structural, copy-number, or representation caveats.** Beyond *CYP2D6*, a small number of panel genes carry more limited caveats:

- ***CYP2B6*** — the common function-defining alleles (e.g., `*6`, `*18`) are SNV-based and are called reliably, but rare *CYP2B6–CYP2B7* hybrid alleles are not resolved by SNV-only calling.
- ***CYP2C19*** — the clinically important alleles (`*2`, `*3`, `*17`, etc.) are SNV-based and fully covered; however, rare copy-number events at the locus are only detectable via CNV analysis (DRAGEN's Star Allele Caller uses CNV calls for *CYP2C19* and *UGT2B17*).
- ***UGT1A1*** — the key alleles `*28` and `*37` (and `*36`) are **promoter TA-repeat length variants** [(TA)₇ and (TA)₈ vs. the (TA)₆ reference; (TA)₅ for `*36`], not simple substitutions. They are callable from a VCF only if the (TA)ₙ indel is correctly represented and left-aligned; mis-representation or low coverage over the TATA box can cause `*28/*37` to be missed. This is a representation/coverage caveat rather than a CNV, but it warrants a targeted check.
- ***G6PD*** — deficiency alleles are SNV-based and are called from the VCF, but *G6PD* is **X-linked**; in hemizygous males a single allele is reported and must be interpreted as such rather than as a homozygous genotype.

Genes not listed here (e.g., *CYP2C9*, *TPMT*, *NUDT15*, *DPYD*, *SLCO1B1*, *VKORC1*, *CYP4F2*, *CYP3A5*, *NAT2*, *F2*, *F5*, *ABCG2*, *IFNL3*, *RYR1*, *CACNA1S*, *CFTR*) are defined by SNVs/indels within the force-called position set and carry no specific SV/CNV limitation in this design.

**6.3 Recommended orthogonal cross-check with DRAGEN CNV, SV, and targeted callers.** Because DRAGEN v4.3 already generates the aligned CRAM in this pipeline, its structural-aware outputs are available at no additional alignment cost and should be consulted **before a structurally-sensitive call is finalized**. In particular:

- **DRAGEN germline CNV and SV callers** can flag copy-number gains/losses and structural breakpoints at the loci above (most importantly the *CYP2D6*/*CYP2D7* region on chromosome 22). A CNV/SV signal there is a direct indication that the SNV-only *CYP2D6* diplotype may be incomplete (e.g., an undetected `*5` deletion or an `xN` duplication) and should be re-adjudicated.
- **DRAGEN targeted callers.** DRAGEN v4.3 additionally ships a **dedicated *CYP2D6* caller derived from the Cyrius method** (determining combined *CYP2D6/CYP2D7* copy number from read depth, detecting SV breakpoints, and emitting structural-aware star-allele genotypes including `xN` duplications and `*5` deletions), a **CYP2B6 caller**, and a **Star Allele Caller** that resolves *CYP2C19* and *UGT2B17* via CNV. Where WGS at ≥30× is available, enabling these callers provides a structurally-complete genotype directly.

The recommended operating procedure is therefore: treat the SNV-only *CYP2D6* diplotype as provisional; consult DRAGEN's CNV/SV output (and, where enabled, the targeted *CYP2D6* caller) for any copy-number or breakpoint signal at the locus; and only then finalize the phenotype.

**6.4 Planned integration path.** The clean, additive route to full structural resolution is to supply a structurally-aware *CYP2D6* diplotype—from DRAGEN's Cyrius-derived targeted caller (or standalone Cyrius; reported concordance ~96.5%, improved to ~99.3% against truth data, versus ~84–87% for prior methods)—to PharmCAT as an **outside call**. PharmCAT natively ingests outside calls, which take precedence over VCF-derived genotypes for that gene, so this augments rather than disrupts the existing workflow. The same mechanism can carry the DRAGEN targeted *CYP2B6* and Star-Allele-Caller (*CYP2C19*, *UGT2B17*) results where structural resolution is desired.

**6.5 Guideline and definition currency.** Interpretive output is only as current as the CPIC guideline set and PharmCAT allele definitions bundled with PharmCAT v3.2.0. Guideline updates, newly defined star alleles, and revised allele-function assignments require a version bump of the interpretation stage; results should always be read together with the pinned tool and definition versions.

**6.6 Analytical, not clinical, scope.** The module reports genotype-derived phenotypes and guideline-concordant annotations. It does not incorporate non-genetic determinants of drug response (co-medication, organ function, adherence, drug–drug interactions), and its output is intended to support—not replace—clinical judgment.

---

## 7. Discussion

The distinguishing engineering choice of the Ilyome PGx module is the elevation of **force-calling to a first-class design principle** rather than a post-hoc patch. Most sequencing-based PGx failures we observe in practice are not errors of the star-allele matcher but errors of *input completeness*: the matcher behaves correctly given ambiguous input, but the input was ambiguous because reference positions were never emitted. By genotyping the entire PGx position set explicitly—upstream of, and complementary to, the PharmCAT preprocessor—the pipeline converts the majority of would-be no-calls into confident reference calls, which raises effective call rates and reduces the number of results a clinician must adjudicate as "inconclusive."

Pairing DRAGEN v4.3 for accelerated, accurately-mapped alignment with GATK HaplotypeCaller v4.6.1 for assembly-based, indel-aware force-calling gives a substrate that is both fast and correct at the positions that matter, while delegating all guideline interpretation to a community-standard, independently-published tool (PharmCAT/CPIC). This separation of concerns—platform-grade alignment, targeted and complete genotyping, and standards-based interpretation—keeps the clinically-sensitive interpretive layer transparent and auditable, and lets the panel track CPIC and PharmVar as they evolve.

The principal remaining gap, *CYP2D6* structural resolution, is well-understood and has a clear, additive solution (Cyrius outside calls). Closing it, together with populating the analytical-validation table against GeT-RM/1000 Genomes reference materials, is the natural next step toward a fully clinically-deployable configuration.

---

## 8. Conclusion

The Ilyome pharmacogenomics module implements a reproducible, guideline-concordant PGx workflow: DRAGEN v4.3 alignment on GRCh38, position-targeted force-calling with GATK HaplotypeCaller v4.6.1 across **1,206 PGx allele-definition positions spanning 22 pharmacogenes**, PharmCAT preprocessing, and PharmCAT v3.2.0 interpretation. By making complete positional genotyping an explicit design goal, the pipeline directly addresses the missing-position problem that undermines call completeness in conventional variant-only workflows, delivering diplotype, allele-function, phenotype, and CPIC-concordant drug guidance in a single auditable report. With the disclosed *CYP2D6* SNV-only limitation and a defined validation and integration roadmap, the module provides a sound foundation for clinical pharmacogenomic reporting on the Ilyome platform.

---

## References

1. Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium of the Pharmacogenomics Research Network. *Clin Pharmacol Ther.* 2011;89(3):464–467. doi:10.1038/clpt.2010.279
2. Relling MV, Klein TE, Gammal RS, Whirl-Carrillo M, Hoffman JM, Caudle KE. The Clinical Pharmacogenetics Implementation Consortium: 10 Years Later. *Clin Pharmacol Ther.* 2020;107(1):171–175. doi:10.1002/cpt.1651
3. Caudle KE, Whirl-Carrillo M, Relling MV, et al. Advancing Clinical Pharmacogenomics Worldwide Through the Clinical Pharmacogenetics Implementation Consortium (CPIC). *Clin Pharmacol Ther.* 2025;118(6):1512–1522. doi:10.1002/cpt.70005
4. Whirl-Carrillo M, McDonagh EM, Hebert JM, et al. Pharmacogenomics Knowledge for Personalized Medicine. *Clin Pharmacol Ther.* 2012;92(4):414–417. doi:10.1038/clpt.2012.96
5. Gong L, Whirl-Carrillo M, Klein TE. PharmGKB, an Integrated Resource of Pharmacogenomic Knowledge. *Curr Protoc.* 2021;1(8):e226. doi:10.1002/cpz1.226
6. Klein TE, Ritchie MD. PharmCAT: A Pharmacogenomics Clinical Annotation Tool. *Clin Pharmacol Ther.* 2018;104(1):19–22. doi:10.1002/cpt.928
7. Sangkuhl K, Whirl-Carrillo M, Whaley RM, et al. Pharmacogenomics Clinical Annotation Tool (PharmCAT). *Clin Pharmacol Ther.* 2020;107(1):203–210. doi:10.1002/cpt.1568
8. Li B, Sangkuhl K, Keat K, et al. How to Run the Pharmacogenomics Clinical Annotation Tool (PharmCAT). *Clin Pharmacol Ther.* 2023;113(5):1036–1047. doi:10.1002/cpt.2790
9. Poplin R, Ruano-Rubio V, DePristo MA, et al. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv.* 2018:201178. doi:10.1101/201178
10. McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res.* 2010;20(9):1297–1303. doi:10.1101/gr.107524.110
11. Behera S, Catreux S, Rossi M, et al. Comprehensive genome analysis and variant detection at scale using DRAGEN. *Nat Biotechnol.* 2024. doi:10.1038/s41587-024-02382-1
12. Gaedigk A, Ingelman-Sundberg M, Miller NA, et al.; PharmVar Steering Committee. The Pharmacogene Variation (PharmVar) Consortium: Incorporation of the Human Cytochrome P450 (CYP) Allele Nomenclature Database. *Clin Pharmacol Ther.* 2018;103(3):399–401. doi:10.1002/cpt.910
13. Chen X, Shen F, Gonzaludo N, et al. Cyrius: accurate CYP2D6 genotyping using whole-genome sequencing data. *Pharmacogenomics J.* 2021;21(2):251–261. doi:10.1038/s41397-020-00205-5
14. Gaedigk A, Turner A, Everts RE, et al. Characterization of Reference Materials for Genetic Testing of CYP2D6 Alleles: A GeT-RM Collaborative Project. *J Mol Diagn.* 2019;21(6):1034–1052. doi:10.1016/j.jmoldx.2019.06.007
15. Gonsalves SG, Dirksen RT, Sangkuhl K, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for the Use of Potent Volatile Anesthetic Agents and Succinylcholine in the Context of RYR1 or CACNA1S Genotypes. *Clin Pharmacol Ther.* 2019;105(6):1338–1344. doi:10.1002/cpt.1319
16. Clancy JP, Johnson SG, Yee SW, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for Ivacaftor Therapy in the Context of CFTR Genotype. *Clin Pharmacol Ther.* 2014;95(6):592–597. doi:10.1038/clpt.2014.54
17. Desta Z, Gammal RS, Gong L, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for CYP2B6 and Efavirenz-Containing Antiretroviral Therapy. *Clin Pharmacol Ther.* 2019;106(4):726–733. doi:10.1002/cpt.1477
18. Moriyama B, Obeng AO, Barbarino J, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for CYP2C19 and Voriconazole Therapy. *Clin Pharmacol Ther.* 2017;102(1):45–51. doi:10.1002/cpt.583
19. Turner AJ, Sangkuhl K, Klein TE, Gaedigk A. PharmVar Tutorial on CYP2D6 Structural Variation Testing and Recommendations on Reporting. *Clin Pharmacol Ther.* 2023;114(5):1032–1043. doi:10.1002/cpt.3044
20. Illumina. DRAGEN v4.3 Product Guide — Targeted Caller (CYP2D6 Caller; Star Allele Caller for CYP2C19 and UGT2B17). Illumina Connected Software documentation, 2024. https://help.dragen.illumina.com

*Additional gene-specific CPIC guidelines relevant to the panel (e.g., DPYD/fluoropyrimidines, TPMT–NUDT15/thiopurines, CYP2C9–VKORC1–CYP4F2/warfarin, SLCO1B1/statins, UGT1A1, G6PD, CYP3A5/tacrolimus, IFNL3) are maintained and versioned by CPIC and are available at cpicpgx.org; the specific guideline version bundled with PharmCAT v3.2.0 governs the reported recommendations.*
