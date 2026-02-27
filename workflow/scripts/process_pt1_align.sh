#!/bin/bash
#SBATCH -t 16:00:00
#SBATCH --mem=35G
#SBATCH --account=hansengroup
#SBATCH -J align
#SBATCH -p himem
#SBATCH -c 6
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS

# Part 1: QC, Trimming, and STAR Alignment
# Submit Part 2 with: sbatch --dependency=afterok:$SLURM_JOB_ID process_part2_quant.sh

set +o pipefail

echo "PWD: $(pwd)"; hostname; date
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"

PIPELINE_STATUS=0

#########################
# Load modules
#########################
module purge
module load fastqc
module load samtools
module load STAR/2.6.1c
module load trim_galore

#########################
# Parameters
#########################
THREADS=10
WORK_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
DB=/cluster/projects/livingbank/workspace/references/hg_38

STAR_INDEX=${DB}/star_index
GTF=${DB}/gencode.v25.annotation.gtf

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
MATE1=$(echo "$SAMPLE_LINE" | cut -d',' -f2)
MATE2=$(echo "$SAMPLE_LINE" | cut -d',' -f3)
GROUP=$(echo "$SAMPLE_LINE" | cut -d',' -f4)

if [[ -z "$NAME" || -z "$MATE1" || -z "$MATE2" ]]; then
    echo "ERROR: Could not parse manifest line $LINE_NUM" >&2
    echo "Line content: $SAMPLE_LINE" >&2
    exit 1
fi

mkdir -p "${WORK_DIR}/${NAME}"
cd "${WORK_DIR}/${NAME}"

echo "============================================"
echo "PART 1: QC, Trimming, and Alignment"
echo "Sample Name: ${NAME}"
echo "Group: ${GROUP}"
echo "Array Task ID: ${IDX}"
echo "Mate 1: ${MATE1}"
echo "Mate 2: ${MATE2}"
echo "Working directory: $(pwd)"
echo "============================================"

rename_log_files() {
    if [[ -n "$SLURM_JOB_ID" && -n "$SLURM_ARRAY_TASK_ID" ]]; then
        LOG_DIR="${WORK_DIR}/logs"
        mkdir -p "$LOG_DIR"
        ORIG_OUT="${LOG_DIR}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
        ORIG_ERR="${LOG_DIR}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
        DESC_OUT="${LOG_DIR}/${NAME}_part1.out"
        DESC_ERR="${LOG_DIR}/${NAME}_part1.err"
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

# Check input files exist
for f in "$MATE1" "$MATE2"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: Missing input: $f" >&2
    exit 1
  fi
done

#########################
# FastQC (raw)
#########################
echo "--- FastQC (raw) ---"
if fastqc -o ./ -t $THREADS "$MATE1" "$MATE2"; then
    echo "FastQC completed successfully"
else
    echo "WARNING: FastQC failed, continuing anyway" >&2
    PIPELINE_STATUS=1
fi

#########################
# Trim Galore (paired, gz)
#########################
echo "--- Trim Galore ---"
if trim_galore --phred33 --paired --output_dir ./ "$MATE1" "$MATE2"; then
    echo "Trim Galore completed successfully"
else
    echo "ERROR: Trim Galore failed" >&2
    exit 1
fi

MATE1_BASE=$(basename "$MATE1")
MATE2_BASE=$(basename "$MATE2")
MATE1_CORE=$(echo "$MATE1_BASE" | sed 's/\.gz$//' | sed 's/\.\(fastq\|fq\)$//')
MATE2_CORE=$(echo "$MATE2_BASE" | sed 's/\.gz$//' | sed 's/\.\(fastq\|fq\)$//')
MATE1_TRIM="./${MATE1_CORE}_val_1.fq.gz"
MATE2_TRIM="./${MATE2_CORE}_val_2.fq.gz"

TRIM_OK=true
for f in "$MATE1_TRIM" "$MATE2_TRIM"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: Trimmed file missing: $f" >&2
    TRIM_OK=false
  fi
done

if [[ "$TRIM_OK" == "false" ]]; then
    echo "ERROR: Cannot continue without trimmed files" >&2
    exit 1
fi

echo "Trimmed reads: $MATE1_TRIM , $MATE2_TRIM"

if zcat "$MATE1_TRIM" | head -n 4 > /dev/null 2>&1; then
    echo "R1 trimmed file validated"
else
    echo "WARNING: R1 trimmed file may be corrupted" >&2
    PIPELINE_STATUS=1
fi

if zcat "$MATE2_TRIM" | head -n 4 > /dev/null 2>&1; then
    echo "R2 trimmed file validated"
else
    echo "WARNING: R2 trimmed file may be corrupted" >&2
    PIPELINE_STATUS=1
fi

#########################
# STAR (genome align + transcriptome BAM + gene counts)
#########################
echo "--- STAR alignment ---"

if STAR \
  --runThreadN $THREADS \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$MATE1_TRIM" "$MATE2_TRIM" \
  --readFilesCommand zcat \
  --sjdbGTFfile "$GTF" \
  --outFileNamePrefix "${NAME}_" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMorder Paired \
  --outBAMsortingThreadN 4 \
  --outBAMsortingBinsN 100 \
  --outSAMattributes All \
  --quantMode TranscriptomeSAM GeneCounts; then
    echo "STAR alignment completed successfully"
else
    echo "ERROR: STAR alignment failed" >&2
    exit 1
fi

BAM_COOR=${NAME}_Aligned.sortedByCoord.out.bam
BAM_TX=${NAME}_Aligned.toTranscriptome.out.bam
LOG_FINAL=${NAME}_Log.final.out

if [[ -s "$BAM_COOR" && -s "$BAM_TX" && -s "$LOG_FINAL" ]]; then
    samtools index -@ "$THREADS" "$BAM_COOR" 2>&1 || echo "WARNING: BAM indexing failed" >&2

    if samtools quickcheck -v "$BAM_COOR" "$BAM_TX" 2>&1; then
        echo "Coordinate and transcriptome BAM validated successfully"
    else
        echo "ERROR: BAM validation failed for sample $NAME" >&2
        exit 1
    fi
else
    echo "ERROR: STAR did not produce required outputs" >&2
    exit 1
fi

#########################
# Strandedness from STAR's ReadsPerGene
#########################
echo "--- Inferring strandedness from STAR ReadsPerGene ---"
RPG=${NAME}_ReadsPerGene.out.tab

if [[ -s "$RPG" ]]; then
    COL2=$(tail -n +5 "$RPG" | awk '{sum+=$2} END{print (sum==""?0:sum)}')
    COL3=$(tail -n +5 "$RPG" | awk '{sum+=$3} END{print (sum==""?0:sum)}')

    STRAND_OPT=no
    if awk "BEGIN{exit !($COL2 > 5*$COL3)}"; then
      STRAND_OPT=reverse
    elif awk "BEGIN{exit !($COL3 > 5*$COL2)}"; then
      STRAND_OPT=yes
    fi
    echo "HTSeq strandedness inferred: ${STRAND_OPT}"
    # Save strandedness for Part 2
    echo "$STRAND_OPT" > "${NAME}_strandedness.txt"
else
    echo "WARNING: Missing ${RPG}, using default strandedness: no" >&2
    STRAND_OPT=no
    echo "$STRAND_OPT" > "${NAME}_strandedness.txt"
    PIPELINE_STATUS=1
fi

#########################
# Mapping summary (partial)
#########################
SUMMARY=${NAME}_mapping_summary.txt
{
  echo "############  Combined R1/R2 mapping summary  #############"
  echo "Sample Name: ${NAME}"
  echo "Group: ${GROUP}"
  echo ""
  if [[ -s "$LOG_FINAL" ]]; then
      cat "$LOG_FINAL"
  else
      echo "STAR Log.final.out not found"
  fi
  echo
  echo "############ Strandedness (from STAR ReadsPerGene) #############"
  echo "HTSeq -s option: ${STRAND_OPT}"
} > "$SUMMARY"

echo "============================================"
echo "PART 1 COMPLETED for sample ${NAME}"
echo "Output BAMs: $BAM_COOR, $BAM_TX"
echo "Run Part 2 for quantification"
echo "============================================"

exit $PIPELINE_STATUS
