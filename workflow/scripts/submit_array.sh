#!/bin/bash

# Script to determine array size and submit the RNA-seq processing job
# This script reads the manifest_combined.csv file which contains:
#   sample_name, mate1, mate2, group
# and submits the array job with the correct parameters

# Work directory where processing happens
WORK_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"

# Default manifest file - can be overridden with first argument
CSV_FILE="${1:-${WORK_DIR}/manifest_combined.csv}"

if [[ ! -f "$CSV_FILE" ]]; then
    echo "Error: CSV file $CSV_FILE not found!" >&2
    echo "Usage: $0 [manifest_combined.csv]" >&2
    exit 1
fi

echo "Analyzing samples in $CSV_FILE..."
echo ""

# Count total samples (excluding header)
total_samples=$(tail -n +2 "$CSV_FILE" | grep -v '^$' | wc -l)

if [[ $total_samples -eq 0 ]]; then
    echo "No samples to process!" >&2
    exit 1
fi

# Count samples by group
echo "Samples by group:"
tail -n +2 "$CSV_FILE" | grep -v '^$' | cut -d',' -f4 | sort | uniq -c | while read count group; do
    echo "  $group: $count samples"
done
echo ""

echo "Summary:"
echo "  Total samples to process: $total_samples"
echo ""

# Generate the sbatch command
array_range="0-$((total_samples-1))"
echo "Suggested SLURM submission command:"
echo "sbatch --array=$array_range process.sh"
echo ""
echo "Alternative: Submit with a limit on concurrent jobs (e.g., max 10 at once):"
echo "sbatch --array=$array_range%10 process.sh"
echo ""

# Show sample information for verification (first 15)
echo "Sample array indices (first 15):"
sample_count=0
tail -n +2 "$CSV_FILE" | grep -v '^$' | while IFS=',' read -r sample_name mate1 mate2 group; do
    echo "  [$sample_count] $sample_name (group: $group)"
    ((sample_count++))

    if [[ $sample_count -ge 15 ]]; then
        break
    fi
done

if [[ $total_samples -gt 15 ]]; then
    echo "  ... (and $((total_samples-15)) more)"
fi
echo ""

# Verify paired files exist (optional check if on cluster)
echo "Note: Before submitting, ensure all FASTQ files exist."
echo "You can verify with: tail -n +2 $CSV_FILE | cut -d',' -f2,3 | tr ',' '\n' | xargs -I {} ls -la {}"
