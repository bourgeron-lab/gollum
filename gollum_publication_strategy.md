# Gollum Publication Strategy

## Option A — Tool / Methods Paper (Recommended)

### Framing

Gollum as a generalizable, lightweight tool for detecting rings on acrocentric chromosomes from standard short-read WGS data. The PMS/r(22) cohort is the discovery and validation dataset, not the main story. Emphasis on: accessibility (GRCh38 input), performance, generalization across acrocentric chromosomes.

### What the paper needs

- **Validation on r(22):** existing 7 cases + new ones. Well-powered, already mostly done.
- **Generalization to other acrocentrics:** 1-2 cases on r(13), r(14), r(15), or r(21). Even a single confirmed case per chromosome is sufficient for a tool paper.
- **Negative controls:** ~220+ non-ring PMS samples + relatives. Report specificity formally.
- **GRCh38 mode benchmark:** show concordance between full-T2T and GRCh38→local-T2T modes. This is a major selling point — most labs have GRCh38 BAMs.
- **Comparison with existing SV callers:** run Manta, DELLY, ClinSV on the same samples and show they miss acrocentric rings. This is the "why does this tool need to exist" argument.
- **Runtime benchmarks:** per-sample runtime, memory usage, scalability.
- **Clean code release:** pip-installable or conda, documentation, test data, Docker image.

### Results structure

1. GRCh38-only mode vs full-T2T mode (concordance, sensitivity)
2. Detection performance on r(22) (sensitivity, specificity, F1)
3. Generalization to other acrocentric chromosomes
4. Benchmark against standard SV callers (Manta, DELLY, ClinSV)
5. Runtime and resource usage
6. (Optional) Mosaicism detection / quantification

### Target journals

| Journal | Impact | Fit | Notes |
|---|---|---|---|
| **Bioinformatics** (Oxford) | ~6 | Excellent | Natural home for bioinformatics tools. Application note (2 pages) or full paper. Reviewers expect clean benchmarks and available code. |
| **GigaScience** | ~7 | Very good | Strong emphasis on open data/code/reproducibility. Would value Docker image, test datasets, Snakemake workflow. |
| **NAR Genomics & Bioinformatics** | ~4 | Good | Accepts focused tool papers. Less competitive, faster review. |
| **Genome Biology** | ~13 | Stretch | Would need broader validation, larger cohorts, or a compelling biological insight beyond the tool itself. |
| **Bioinformatics Advances** | ~3 | Good | Newer OUP journal, faster turnaround, lower bar than Bioinformatics. |

### Pros

- Aligns with the planned heavy tool update (minimap2, GRCh38 mode, HDBSCAN, etc.)
- Limited cohort size is less of an issue — tool papers are judged on methodology and benchmarks
- High reuse value — others can apply it to their own cohorts
- Clean narrative: "existing pipelines miss acrocentric rings, here's a tool that finds them"

### Cons

- Need to invest in engineering (packaging, documentation, Docker)
- Benchmark against other tools requires running multiple pipelines
- Reviewers may ask for long-read comparison or simulation studies

---

## Option B — Clinical Genomics Paper

### Framing

WGS combined with T2T reference reveals underdiagnosed ring chromosomes in PMS patients, with direct clinical implications for NF2 surveillance. Gollum is presented as the method enabling the finding, not as the main contribution.

### What the paper needs

- **Expanded PMS cohort with clinical data:** detailed phenotyping, NF2 screening results, longitudinal follow-up.
- **New r(22) diagnoses:** the P7 discovery story is compelling — more cases like this strengthen the clinical message.
- **NF2 risk quantification:** ideally, data on tumor incidence in r(22) vs simple deletion PMS patients. Literature review if own data is insufficient.
- **Genotype-phenotype correlations:** ring vs non-ring PMS patients, severity scores, speech outcomes, etc.
- **Gollum as supplementary method:** described in methods/supplementary, not the focus.

### Results structure

1. Cohort description and genetic findings
2. Ring chromosome detection using WGS + gollum
3. Complex ring characterization (P1 inv-dup-del)
4. Clinical comparison: ring vs non-ring PMS patients
5. NF2 risk and surveillance recommendations

### Target journals

| Journal | Impact | Fit | Notes |
|---|---|---|---|
| **npj Genomic Medicine** | ~7 | Good | Current target. Reasonable fit for genomics + clinical implications. |
| **Genetics in Medicine** (ACMG) | ~9 | Very good | Strong clinical genetics audience. Would want clear clinical actionability (NF2 surveillance protocol). |
| **European J Human Genetics** | ~5 | Good | European cohort, clinical focus, solid readership in rare disease genetics. |
| **AJHG** | ~11 | Stretch | Would need larger cohort, multiple chromosome types, and a strong mechanistic story. |
| **Human Mutation / Human Genomics** | ~3-4 | Fallback | Solid mid-tier, less competitive. |

### Pros

- More directly clinically impactful (NF2 surveillance recommendations)
- Leverages the existing clinical data and phenotyping work
- The P1 complex ring is a compelling case study

### Cons

- Small cohort (7-12 ring cases) limits statistical power for genotype-phenotype analysis
- Gollum improvements are undervalued — the tool is just a means to an end
- Harder to differentiate from existing PMS literature without a larger cohort
- Clinical reviewers may want longer follow-up data or prospective NF2 screening

---

## Recommendation

**Go with Option A (tool paper)** given:

1. You're planning substantial tool updates — this work deserves to be the main contribution
2. Cohort size (10-12 cases) is sufficient for a tool validation but weak for clinical genomics
3. The GRCh38 mode is a concrete, publishable innovation that doesn't exist elsewhere
4. A separate clinical brief communication / letter can follow later with more cases and follow-up

**Fallback plan:** if the tool updates take longer than expected or validation cases don't materialize, Option B with the existing data + updated gollum in supplementary is a viable path to publication with less new work required.
