# Gollum v2 — Technical Implementation Roadmap

## 1. Core Pipeline Redesign

### 1.1 Replace BLAT with minimap2

**Current:** BLAT is used to realign R2 reads against T2T-CHM13. Slow, aging, licensing issues for redistribution.

**New:** minimap2 with short-read preset.

```bash
minimap2 -a -x sr --secondary=yes -N 5 t2t_chm13.fa r2_reads.fq
```

- `-x sr` for short reads (150bp Illumina)
- `--secondary=yes -N 5` to report multi-mapping (critical for PHR regions)
- Output SAM for direct access to AS, NM, de tags

**Filtering strategy:**

minimap2 does not output e-values or bit scores. Replace with:

| Metric | Source | Threshold | Purpose |
|---|---|---|---|
| `AS` (alignment score) | SAM tag | ≥ 100 | Replaces BLAT bit score. Scales with read length. |
| `de` (sequence divergence) | SAM tag | ≤ 0.05 | Replaces e-value. Max 5% divergence. |
| Location | SAM RNAME/POS | Within any SAAC region | Same as current approach. |

**Important:** Do NOT filter R2 on `mapQ`. Reads mapping to PHRs will have mapQ=0 due to multi-mapping — this is expected and correct. Breakpoint specificity comes from the R1 anchor (which maps uniquely to the q-arm), not from R2 disambiguation.

**Optional — specificity ratio for non-PHR reads:**

```python
# From minimap2 SAM tags
# s1 = best chain score, s2 = second-best chain score
if s2 / s1 < 0.8:  # read maps specifically to one SAAC
    # Informative: can narrow down which short arm is involved
    # Mostly useful outside PHR regions
```

In practice, most R2 reads will be in PHRs (ratio ≈ 1.0), so this is secondary information.

**Validation:** Run both BLAT and minimap2 on existing r(22) cases. Compare number of R2 reads passing filters. Expect similar or better sensitivity with minimap2.

---

### 1.2 Clustering: HDBSCAN single-pass (recommended) vs two-pass DBSCAN

#### Current approach (two-pass DBSCAN)

1. DBSCAN on raw R1 positions → candidate regions
2. Extract R2, realign with BLAT, filter
3. DBSCAN again on filtered R1 positions → final breakpoints

**Problems:**
- First DBSCAN may discard low-signal cases (mosaic rings, low coverage) before realignment
- Same eps/min_samples used in different noise contexts
- No confidence score — binary cluster membership

#### Option A — Single-pass after minimap2 (recommended)

With minimap2 replacing BLAT, realignment is fast enough to skip pre-filtering:

1. Extract ALL discordant reads (R1 on target q-arm, R2 on any SAAC or unmapped)
2. Realign ALL R2 with minimap2
3. Filter on AS/de + SAAC location
4. Single HDBSCAN on surviving R1 positions

**Advantages:** simpler, no risk of losing low-signal cases, fewer parameters.

**Runtime cost:** negligible — typically hundreds to low thousands of candidate R2 reads per sample, minimap2 handles this in seconds.

#### Option B — Keep two-pass but improve

If keeping two passes for any reason:
- Use HDBSCAN instead of DBSCAN for both passes
- Set `min_samples` adaptively based on coverage: `min_samples = max(3, int(coverage / 15))`
- Use different parameters for each pass (first pass: permissive, second pass: stricter)

#### HDBSCAN vs DBSCAN

| Feature | DBSCAN | HDBSCAN |
|---|---|---|
| Requires eps parameter | Yes | No |
| Handles variable density | Poorly | Well |
| Confidence scores | No | Yes (probabilities per point) |
| Noise detection | Binary | Soft (outlier scores) |
| Install | sklearn built-in | `pip install hdbscan` |

```python
from hdbscan import HDBSCAN

clusterer = HDBSCAN(min_cluster_size=3, min_samples=2)
labels = clusterer.fit_predict(r1_positions.reshape(-1, 1))
probabilities = clusterer.probabilities_  # per-read confidence

# Report per-cluster:
# - number of supporting reads
# - mean probability
# - genomic span
```

**Recommendation:** Single-pass + HDBSCAN. Report cluster size + mean probability as confidence score per candidate breakpoint.

---

### 1.3 GRCh38 Input Mode

**Rationale:** Most clinical labs have GRCh38-aligned BAM/CRAM. Requiring full T2T realignment is the biggest adoption barrier.

**Workflow:**

```
GRCh38 BAM/CRAM
    │
    ├─ Extract discordant reads from target q-arm (samtools view)
    │   - R2 unmapped OR mapped to chrUn/random/alt contigs
    │   - R2 with low mapQ (multi-mapping)
    │
    ├─ Also extract reads near the terminal deletion breakpoint
    │   - Soft-clipped reads at the breakpoint
    │   - Reads with unmapped mates in the breakpoint region
    │
    └─ Realign extracted R2 to T2T-CHM13 with minimap2
        │
        └─ Standard filtering + HDBSCAN (same as full-T2T mode)
```

**Key difference from full-T2T mode:** In GRCh38, SAAC sequences are absent, so R2 reads from ring breakpoints will be:
- Unmapped
- Mapped to `chrUn_*` or alt contigs with low quality
- Mapped elsewhere with poor alignment

The extraction step needs to capture all of these. Use samtools flags:

```bash
# Extract R1 on target q-arm with R2 unmapped or poorly mapped
samtools view -h input.cram target_region | \
    awk '$7 == "*" || $5 < 10 || $7 ~ /chrUn/' > candidates.sam
```

**Validation needed:** Run GRCh38 mode and full-T2T mode on same samples. Report:
- Concordance (same breakpoints detected?)
- Sensitivity difference (any missed calls in GRCh38 mode?)
- Runtime comparison

**This is a publishable result on its own** — demonstrating that acrocentric ring detection is possible from standard GRCh38 data without full genome realignment.

---

## 2. Scope Extension

### 2.1 All Acrocentric Chromosomes (13, 14, 15, 21, 22)

**Current:** Designed and tested on chr22 only.

**Changes needed:**
- Parameterize the target chromosome (already partially done?)
- Define SAAC regions for each acrocentric on T2T-CHM13
- Define q-arm boundaries for each chromosome
- Test on available cases (1-2 per chromosome)

**SAAC region definitions on T2T-CHM13:**

These need to be extracted from the T2T annotation. For each acrocentric, the short arm spans from position 0 to the centromere boundary. The PHR regions (Guarracino et al. 2023) overlap across all 5 chromosomes.

```python
ACROCENTRIC_CONFIG = {
    "chr13": {"q_start": ..., "q_end": ..., "p_start": 0, "p_end": ...},
    "chr14": {"q_start": ..., "q_end": ..., "p_start": 0, "p_end": ...},
    "chr15": {"q_start": ..., "q_end": ..., "p_start": 0, "p_end": ...},
    "chr21": {"q_start": ..., "q_end": ..., "p_start": 0, "p_end": ...},
    "chr22": {"q_start": ..., "q_end": ..., "p_start": 0, "p_end": ...},
}
```

**Complication:** Each chromosome may have different noise profiles and blacklist regions. The chr22 blacklist won't transfer directly.

### 2.2 Non-acrocentric Chromosomes (exploratory)

If you get 1-2 cases on non-acrocentric chromosomes, the detection logic is fundamentally different:
- Both breakpoints are on mappable arms
- Standard SV callers (Manta, DELLY) *should* detect these
- Gollum's value-add is unclear for non-acrocentrics

**Recommendation:** Test standard SV callers on these cases. If they detect the ring → mention in discussion as context. If they miss it → interesting finding, but likely edge case.

Don't build non-acrocentric support into gollum unless there's a clear gap that existing tools miss.

---

## 3. Benchmarking

### 3.1 Against Existing SV Callers

**Goal:** Demonstrate that standard SV callers fail to detect acrocentric rings.

**Tools to test:**

| Tool | Type | Why include |
|---|---|---|
| Manta | SV caller (paired-end + split-read) | Most widely used clinical SV caller |
| DELLY | SV caller (paired-end + split-read) | Already used in your pipeline |
| ClinSV | Clinical SV caller | Already used in your pipeline, clinical-grade |
| Smoove | SV caller (LUMPY wrapper) | Common in research pipelines |

**Run each on:**
- All r(22) samples (positive controls)
- Non-ring PMS deletion samples (negative controls)
- Both GRCh38 and T2T-CHM13 alignments

**Metrics to report:**

| Metric | Definition |
|---|---|
| Detection rate | Does the tool output any call overlapping the ring breakpoint? |
| Call type | What does it call it? (deletion, translocation, BND, nothing?) |
| Breakpoint accuracy | Distance from called breakpoint to true breakpoint (bp) |

**Expected result:** All tools detect the terminal deletion. None detect the ring (junction to short arm) because GRCh38 lacks SAAC sequences. On T2T, some may partially detect it as a translocation/BND, but without the ring interpretation.

### 3.2 Gollum Internal Benchmarks

| Metric | How to compute |
|---|---|
| Sensitivity | TP / (TP + FN) across all confirmed ring cases |
| Specificity | TN / (TN + FP) across all non-ring controls |
| F1 score | Harmonic mean of precision and recall |
| Runtime | Wall-clock time per sample (mean ± sd) |
| Memory | Peak RSS per sample |
| GRCh38 vs T2T concordance | % of calls identical between modes |

### 3.3 Simulated Data (optional but recommended)

Simulate ring-supporting read pairs at known positions with varying:
- Coverage (5x, 10x, 15x, 30x, 50x)
- Mosaicism level (5%, 10%, 25%, 50%, 100%)
- Breakpoint location (within PHR, outside PHR, near blacklist)

This gives you sensitivity curves and defines the detection limit. Tools like `wgsim` or `art_illumina` can generate reads; you'd need to engineer the discordant pairs manually.

---

## 4. Additional Features

### 4.1 Quantitative Mosaicism Estimation

**Current:** Binary call (ring / no ring).

**Proposed:** Estimate the fraction of cells carrying the ring.

```python
# Simplified approach
def estimate_mosaicism(n_supporting_reads, coverage_at_breakpoint, read_length=150, insert_size=400):
    """
    Expected supporting reads for 100% ring:
    ~ coverage * (2 * read_length) / insert_size at the breakpoint
    (rough approximation — the actual expectation depends on
    fragment size distribution and breakpoint geometry)
    """
    expected_full = coverage_at_breakpoint * (2 * read_length) / insert_size
    vaf = n_supporting_reads / expected_full
    return min(vaf, 1.0)  # cap at 1.0
```

This is a rough estimate — proper mosaicism quantification would need a statistical model accounting for fragment size distribution, mapping biases near breakpoints, and sampling variance. But even a rough number (e.g., "~50% mosaic") is clinically useful.

**Validation:** Compare with cytogenetic mosaicism estimates (karyotype on X metaphases) if available.

### 4.2 Blacklist / Panel of Normals

**Current:** Ad hoc exclusion of "regions having candidate ring breakpoints across multiple individuals."

**Proposed:** Build a formal panel of normals (PoN).

1. Run gollum on N control samples (non-ring, ideally diverse ancestry)
2. Record all candidate breakpoint positions
3. Any position with hits in ≥ K samples (e.g., K=3) is blacklisted
4. Store as BED file, ship with gollum

```
# blacklist.bed
chr22   16000000   16050000   pon_hit_count=15
chr22   18500000   18520000   pon_hit_count=8
...
```

Use the ~220 non-ring samples from the PMS cohort as initial PoN. Users can extend with their own controls.

### 4.3 Output Format

**Current:** unclear (custom format?).

**Proposed:** Standard BED + optional VCF.

```
# gollum_output.bed
#chrom  start       end         name            score   info
chr22   47006000    47007000    ring_bp_1       25      n_reads=25;confidence=0.95;mosaicism=0.85;mode=t2t
```

VCF with BND notation for interoperability with other SV tools:

```
#CHROM  POS         ID          REF ALT                         FILTER  INFO
chr22   47006904    gollum_1    N   N]chr22p:8500000]           PASS    SVTYPE=BND;SUPPORT=25;CONF=0.95
```

---

## 5. Engineering

### 5.1 Packaging

- pip-installable (`pip install gollum`)
- Conda recipe (bioconda channel)
- Docker/Singularity image with minimap2 + samtools bundled
- Snakemake or Nextflow wrapper for batch processing

### 5.2 Dependencies (new)

| Dependency | Purpose | Replaces |
|---|---|---|
| minimap2 | R2 realignment | BLAT |
| hdbscan | Clustering | sklearn DBSCAN |
| samtools | Read extraction | samtools (unchanged) |
| pysam | BAM/CRAM parsing | — |

### 5.3 Test Suite

- Unit tests for filtering, clustering, mosaicism estimation
- Integration tests with bundled mini-BAM files from known r(22) cases
- CI/CD with GitHub Actions

---

## 6. Priority / Difficulty Matrix

| Change | Difficulty | Priority for publication | Notes |
|---|---|---|---|
| **Replace BLAT with minimap2** | Low | **Critical** | Straightforward swap. Removes licensing issue. Faster. |
| **Single-pass clustering** | Low | **High** | Simplifies pipeline. Direct consequence of minimap2 speed. |
| **HDBSCAN instead of DBSCAN** | Low | **High** | Drop-in replacement. Gives confidence scores for free. |
| **GRCh38 input mode** | Medium | **Critical** | Key selling point. Needs careful extraction logic for unmapped/poorly-mapped R2. Requires validation against full-T2T mode. |
| **Extend to all acrocentrics** | Medium | **Critical** | Parameterize target chromosome + define SAAC regions. Need 1-2 validation cases per chromosome. |
| **Benchmark vs other SV callers** | Medium | **Critical** | Running Manta/DELLY/ClinSV on all samples. Mostly execution time, not complexity. |
| **Blacklist / Panel of normals** | Low-Medium | **High** | Run gollum on controls, aggregate. Replaces ad hoc exclusion. |
| **Output format (BED + VCF)** | Low | **High** | Standardize output. Enables downstream integration. |
| **Quantitative mosaicism** | Medium | **Medium** | Nice to have. Rough estimate is easy; proper model is harder. Clinically relevant but not strictly needed for tool paper. |
| **Simulated data benchmarks** | Medium | **Medium** | Sensitivity curves by coverage/mosaicism. Strengthens the paper but time-consuming to set up properly. |
| **Packaging (pip/conda/Docker)** | Low-Medium | **High** | Expected by reviewers for Bioinformatics/GigaScience. Docker is near-mandatory. |
| **Snakemake/Nextflow wrapper** | Low | **Low** | Nice to have. Not expected for initial publication. |
| **Non-acrocentric support** | Medium-High | **Low** | Only if cases show existing tools fail. Otherwise just discuss. |
| **Test suite + CI/CD** | Low | **Medium** | Good practice. Some journals (GigaScience) value this explicitly. |

### Suggested implementation order

1. **minimap2 swap + single-pass HDBSCAN** (foundation — everything else depends on this)
2. **Extend to all acrocentrics** (parameterize config)
3. **GRCh38 mode** (biggest user-facing improvement)
4. **Blacklist from panel of normals** (run on controls)
5. **Benchmark against SV callers** (run pipelines, collect results)
6. **Output format standardization** (BED + VCF)
7. **Packaging + Docker** (before submission)
8. **Mosaicism estimation** (if time permits)
9. **Simulated data** (if time permits)
10. **Test suite** (ongoing, alongside development)
