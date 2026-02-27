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

### 5. Build expression matrices

Aggregates RSEM output into gene x sample matrices.

```bash
bash combine_rsem.sh
```

Outputs:
- `combined_TPM.txt` — TPM values (recomputed from expected counts and effective length)
- `combined_expected_counts.txt` — use as input for DESeq2/edgeR

### 6. Gene detection QC

```bash
bash combine_gene_counts.sh
```

Output: `gene_detection_stats.tsv` — per-sample count of genes detected at TPM > 0, 0.1, 1, and 10.

### 7. Alignment QC

```bash
bash combine_star.sh
```

Aggregates STAR alignment statistics across samples.

---

## MHC expression

```bash
bash 00_extract_MHC.sh   # extract MHC gene expression
bash 01_rank_MHC.sh      # rank samples by MHC expression
```
