#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --mem=20G
#SBATCH --account=hansengroup
#SBATCH -J quant
#SBATCH -p himem
#SBATCH -c 6
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS

# Part 2: Quantification (HTSeq, RSEM, StringTie)
# Requires Part 1 to have completed successfully
# Submit with: sbatch --dependency=afterok:<PART1_JOB_ID> --array=0-N process_part2_quant.sh

set +o pipefail

echo "PWD: $(pwd)"; hostname; date
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"

PIPELINE_STATUS=0

#########################
# Load modules
#########################
module purge
module load samtools
module load rsem
module load HTSeq/0.11.0
PATH=$PATH:/cluster/home/t117036uhn/local_soft/stringtie2

#########################
# Parameters
#########################
THREADS=10
WORK_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
DB=/cluster/projects/livingbank/workspace/references/hg_38

FA=${DB}/GRCh38.primary_assembly.genome.fa
GTF=${DB}/gencode.v25.annotation.gtf
RSEM_REF_DIR=${DB}/rsem_reference

MANIFEST_FILE="${WORK_DIR}/manifest_combined.csv"

if [[ ! -s "$MANIFEST_FILE" ]]; then
    echo "ERROR: Manifest file not found: $MANIFEST_FILE" >&2
    exit 1
fi

IDX=${SLURM_ARRAY_TASK_ID:?Set --array or export SLURM_ARRAY_TASK_ID}
LINE_NUM=$((IDX + 2))

TOTAL_SAMPLES=$(tail -n +2 "$MANIFEST_FILE" | grep -v '^$' | wc -l)

if [[ $IDX -ge $TOTAL_SAMPLES ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($IDX) exceeds number of samples ($TOTAL_SAMPLES)" >&2
    exit 1
fi

SAMPLE_LINE=$(sed -n "${LINE_NUM}p" "$MANIFEST_FILE")
NAME=$(echo "$SAMPLE_LINE" | cut -d',' -f1)
GROUP=$(echo "$SAMPLE_LINE" | cut -d',' -f4)

if [[ -z "$NAME" ]]; then
    echo "ERROR: Could not parse manifest line $LINE_NUM" >&2
    exit 1
fi

cd "${WORK_DIR}/${NAME}"

echo "============================================"
echo "PART 2: Quantification"
echo "Sample Name: ${NAME}"
echo "Group: ${GROUP}"
echo "Array Task ID: ${IDX}"
echo "Working directory: $(pwd)"
echo "============================================"

rename_log_files() {
    if [[ -n "$SLURM_JOB_ID" && -n "$SLURM_ARRAY_TASK_ID" ]]; then
        LOG_DIR="${WORK_DIR}/logs"
        mkdir -p "$LOG_DIR"
        ORIG_OUT="${LOG_DIR}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
        ORIG_ERR="${LOG_DIR}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
        DESC_OUT="${LOG_DIR}/${NAME}_part2.out"
        DESC_ERR="${LOG_DIR}/${NAME}_part2.err"
        if [[ -f "$ORIG_OUT" ]]; then
            cp "$ORIG_OUT" "$DESC_OUT"
            echo "Log copied to: $DESC_OUT"
        fi
        if [[ -f "$ORIG_ERR" ]]; then
            cp "$ORIG_ERR" "$DESC_ERR"
            echo "Error log copied to: $DESC_ERR"
        fi
    fi
}

trap rename_log_files EXIT

# Check that Part 1 outputs exist
BAM_COOR=${NAME}_Aligned.sortedByCoord.out.bam
BAM_TX=${NAME}_Aligned.toTranscriptome.out.bam

if [[ ! -s "$BAM_COOR" || ! -s "$BAM_TX" ]]; then
    echo "ERROR: Part 1 outputs not found. Run process_part1_align.sh first." >&2
    echo "Missing: $BAM_COOR or $BAM_TX" >&2
    exit 1
fi

# Load strandedness from Part 1
STRAND_FILE="${NAME}_strandedness.txt"
if [[ -s "$STRAND_FILE" ]]; then
    STRAND_OPT=$(cat "$STRAND_FILE")
    echo "Loaded strandedness from Part 1: ${STRAND_OPT}"
else
    echo "WARNING: Strandedness file not found, using default: no" >&2
    STRAND_OPT=no
fi

#########################
# HTSeq (coding genes only)
#########################
echo "--- HTSeq-count (coding only) ---"

GTF_CODING=${DB}/gencode.v25.annotation.coding_only.gtf
if [[ ! -s "$GTF_CODING" ]]; then
  echo "Creating coding-only GTF at ${GTF_CODING}"
  if grep 'gene_type "protein_coding"' "$GTF" > "$GTF_CODING"; then
      echo "Coding GTF created successfully"
  else
      echo "WARNING: Failed to create coding GTF" >&2
      PIPELINE_STATUS=1
  fi
fi

if [[ -s "$BAM_COOR" && -s "$GTF_CODING" ]]; then
    if htseq-count \
      -f bam -r pos -s "${STRAND_OPT}" \
      -t exon -i gene_id \
      "$BAM_COOR" "$GTF_CODING" > ${NAME}_count_coding.txt; then
        echo "HTSeq count completed successfully"
    else
        echo "WARNING: HTSeq count failed" >&2
        PIPELINE_STATUS=1
    fi
else
    echo "WARNING: Skipping HTSeq (missing BAM or GTF)" >&2
    PIPELINE_STATUS=1
fi

#########################
# RSEM (from transcriptome BAM)
#########################
echo "--- RSEM quantification ---"

if [[ ! -d "$RSEM_REF_DIR" || ! -s ${RSEM_REF_DIR}/GRCh38.grp ]]; then
  echo "Preparing RSEM reference in ${RSEM_REF_DIR}"
  mkdir -p "$RSEM_REF_DIR"
  if rsem-prepare-reference --gtf "$GTF" --num-threads $THREADS "$FA" ${RSEM_REF_DIR}/GRCh38; then
      echo "RSEM reference prepared successfully"
  else
      echo "WARNING: RSEM reference preparation failed" >&2
      PIPELINE_STATUS=1
  fi
fi

if [[ -s "$BAM_TX" ]]; then
    if rsem-calculate-expression \
      --paired-end \
      --bam \
      --no-bam-output \
      --num-threads $THREADS \
      "$BAM_TX" \
      ${RSEM_REF_DIR}/GRCh38 \
      ${NAME}_rsem; then
        echo "RSEM quantification completed successfully"
    else
        echo "WARNING: RSEM quantification failed" >&2
        PIPELINE_STATUS=1
    fi
else
    echo "WARNING: Skipping RSEM (missing transcriptome BAM)" >&2
    PIPELINE_STATUS=1
fi

#########################
# StringTie (optional)
#########################
echo "--- StringTie (optional) ---"
if command -v stringtie >/dev/null 2>&1; then
  if [[ -s "$BAM_COOR" ]]; then
      if stringtie \
        -o ${NAME}_annotation.gtf \
        -G "$GTF" \
        -v -A ${NAME}_abundance.tab \
        -p $THREADS -B \
        "$BAM_COOR"; then
          echo "StringTie quantification completed successfully"
      else
          echo "WARNING: StringTie failed" >&2
          PIPELINE_STATUS=1
      fi
  else
      echo "WARNING: Skipping StringTie (missing BAM)" >&2
  fi
else
  echo "StringTie not found on PATH; skipping."
fi

#########################
# Summary: counts & top expressed
#########################
SUMMARY=${NAME}_mapping_summary.txt
{
  echo
  echo "############ Gene quantification summary #############"
  if [[ -s ${NAME}_ReadsPerGene.out.tab ]]; then
    TOTAL_GENES=$(tail -n +5 ${NAME}_ReadsPerGene.out.tab | wc -l)
    echo "STAR quantified genes (reads-per-gene table): ${TOTAL_GENES}"
  fi

  if [[ -s ${NAME}_rsem.genes.results ]]; then
    TOTAL_RSEM=$(tail -n +2 ${NAME}_rsem.genes.results | wc -l)
    echo "RSEM genes: ${TOTAL_RSEM}"
    echo "Top 10 genes by TPM:"
    tail -n +2 ${NAME}_rsem.genes.results | sort -k6,6nr | head -10 | cut -f1,6
  fi

  echo
  echo "############ Output file sizes #############"
  shopt -s nullglob
  files=( ${NAME}_*.bam ${NAME}_*.bai ${NAME}_*.gtf ${NAME}_*.tab ${NAME}_*.txt ${NAME}_*.results )
  if ((${#files[@]})); then
    for f in "${files[@]}"; do
      size=$(du -h "$f" | awk '{print $1}')
      echo "$size $f"
    done
  else
    echo "No matching output files found for ${NAME}_*"
  fi
  shopt -u nullglob

  echo
  echo "############ Pipeline Status #############"
  if [[ $PIPELINE_STATUS -eq 0 ]]; then
      echo "All processing steps completed successfully."
  else
      echo "Pipeline completed with WARNINGS/ERRORS. Check log files for details."
  fi
} >> "$SUMMARY"

echo "============================================"
echo "PART 2 COMPLETED for sample ${NAME} (group: ${GROUP})"
echo "============================================"

exit $PIPELINE_STATUS
