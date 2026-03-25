# PMLB Organoid RNA-seq Analysis

RNA-seq processing pipeline for organoid samples from the [Princess Margaret Living Biobank (PMLB)](https://pmlivingbiobank.uhnresearch.ca/). Covers FASTQ concatenation, alignment (STAR), quantification (RSEM), expression matrix generation, QC, and MHC expression extraction.

---

## Data files

| File | Description |
|------|-------------|
| `260225_manifest_replicates.csv` | Per-replicate sample manifest with `group_name` column for concatenation |
| `sample_manifest.tsv` | Complete list of sample replicates |

---

## Pipeline

### 1. Verify source files are accessible from the cluster

```bash
ls /cluster/projects/livingbank/Data/RNAseq/Cescon_lab/ | head -3
ls /cluster/projects/livingbank/Data/RNAseq/RNAseq_Niki_Novogene/ | head -3
ls /cluster/projects/livingbank/Data/RNAseq/PMLB/ | head -3
```

### 2. Concatenate FASTQ files

Groups samples by `group_name`; samples with an empty `group_name` remain separate.

```bash
sbatch concatenate.sh 260225_manifest_replicates.csv manifest_final.tsv
```

- `.out` log prints `Samples written: N` and `Samples skipped: N` at the end
- `.err` log: `grep WARNING logs/<jobid>.out` shows any unreadable source files
- Verify output: `wc -l manifest_final.tsv` should equal expected sample count + 1 (header)

### 3. Generate array submission command

```bash
bash submit_array.sh
# Prints the sbatch --array=0-N command to use in the next step
```

### 4. Alignment + quantification

Submit alignment (Part 1) and quantification (Part 2) as dependent jobs:

```bash
JOB1=$(sbatch --array=0-231 process_pt1_align.sh | awk '{print $4}')
sbatch --dependency=afterok:$JOB1 --array=0-231 process_pt2_quant.sh
```

### 5. Build STAR count matrix

Combines per-sample `ReadsPerGene.out.tab` files into a single gene × sample integer count matrix, selecting the correct strand column per sample based on the strandedness inferred during alignment.

```bash
bash combine_star_counts.sh
```

Output: `combined_star_counts.txt` — rows are ENSG gene IDs, columns are samples, values are integer read counts.

### 6. Alignment QC

```bash
bash combine_star.sh
```

Aggregates STAR alignment statistics across samples into `STAR_mapping_stats.tsv`.

### 7. Batch correction & TPM

Applied in `workflow/analysis/temp_02_batch_correction.ipynb`. See `workflow/analysis/README_batchcorr.txt` for full details.

Samples were first filtered for QC: uniquely mapped % and library size were each compared against the cohort median, and samples falling more than 3 MADs from the median were removed (40 dropped, 192 retained). Genes were then filtered to those with raw count > 10 in at least 5 samples.

Batch correction was performed using ComBat-seq (`inmoose.pycombat`), which models count data with a negative binomial distribution and returns adjusted integer counts. `cancer_type` was included as a biological covariate to prevent over-correction of biological signal. The correction was run independently on two gene sets: all expressed genes (30,552) and protein-coding genes only (17,273), identified from GENCODE v25. Samples belonging to singleton batches were excluded from the model and re-appended uncorrected.

TPM was computed from the protein-coding batch-corrected counts using the longest annotated CDS per gene from the GENCODE v25 GTF (rather than RSEM per-sample effective lengths, which vary across samples). For each gene *g* and sample *s*:

```
RPK(g, s)  = corrected_count(g, s) / (longest_CDS(g) / 1000)
TPM(g, s)  = RPK(g, s) / Σ_g RPK(g, s) × 1,000,000
```

Outputs (`results/`):
- `batch_corrected_counts_all_genes.csv` — adjusted integer counts, 30,552 expressed genes × 192 samples
- `batch_corrected_counts_protein_coding.csv` — adjusted integer counts, 17,273 protein-coding genes × 192 samples
- `batch_corrected_TPM_protein_coding.csv` — TPM from protein-coding batch-corrected counts using longest CDS, 17,273 genes × 192 samples
- `batch_corrected_TPM_all_genes.csv` — TPM from all-genes batch-corrected counts using longest CDS; genes lacking any CDS annotation in GENCODE v25 (lncRNAs, pseudogenes, etc.) are excluded

---

## MHC expression

```bash
bash 00_extract_MHC.sh   # extract MHC gene expression
bash 01_rank_MHC.sh      # rank samples by MHC expression
```
