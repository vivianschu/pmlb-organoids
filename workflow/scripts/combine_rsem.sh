#!/bin/bash
# usage: bash combine_rsem.sh

# Directory containing all sample folders
BASE_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process" 
OUTPUT_DIR="${BASE_DIR}/combined_results"

mkdir -p "$OUTPUT_DIR"

# Find all RSEM gene results files
# Files are in ${BASE_DIR}/{sample_name}/{sample_name}_rsem.genes.results
RSEM_FILES=($(find "$BASE_DIR" -maxdepth 2 -name "*_rsem.genes.results" -type f | sort))

if [ ${#RSEM_FILES[@]} -eq 0 ]; then
    echo "No RSEM gene results files found!"
    exit 1
fi

echo "Found ${#RSEM_FILES[@]} samples"

# Create TPM matrix
# TPM is computed from expected_count (col 5) and effective_length (col 4):
#   RPK = expected_count / (effective_length / 1000)
#   TPM = (RPK / sum(RPK)) * 1e6
echo "Creating TPM matrix..."
# Start with gene_id column from first file
cut -f1 "${RSEM_FILES[0]}" > "$OUTPUT_DIR/gene_ids.tmp"

for file in "${RSEM_FILES[@]}"; do
    sample=$(basename "$file" "_rsem.genes.results")
    echo "$sample"
    awk -F'\t' -v sample="$sample" '
        NR==1 { print sample; next }
        {
            count[NR] = $5
            eff_len[NR] = $4
            rpk[NR] = (eff_len[NR] > 0) ? count[NR] / (eff_len[NR] / 1000) : 0
            sum_rpk += rpk[NR]
        }
        END {
            for (i=2; i<=NR; i++) {
                tpm = (sum_rpk > 0) ? (rpk[i] / sum_rpk) * 1e6 : 0
                printf "%.6f\n", tpm
            }
        }
    ' "$file" > "$OUTPUT_DIR/${sample}_tpm.tmp"
done

# Paste all together
paste "$OUTPUT_DIR/gene_ids.tmp" "$OUTPUT_DIR"/*_tpm.tmp > "$OUTPUT_DIR/combined_TPM.txt"

# Create expected counts matrix (useful for DESeq2)
echo "Creating expected counts matrix..."
for file in "${RSEM_FILES[@]}"; do
    sample=$(basename "$file" "_rsem.genes.results")
    cut -f5 "$file" | sed "1s/expected_count/$sample/" > "$OUTPUT_DIR/${sample}_counts.tmp"
done

paste "$OUTPUT_DIR/gene_ids.tmp" "$OUTPUT_DIR"/*_counts.tmp > "$OUTPUT_DIR/combined_expected_counts.txt"

# Cleanup temp files
rm -f "$OUTPUT_DIR"/*.tmp

echo "Done! Output files:"
echo "  - $OUTPUT_DIR/combined_TPM.txt"
echo "  - $OUTPUT_DIR/combined_expected_counts.txt"

# Show preview
echo ""
echo "Preview of TPM matrix (first 5 lines, first 5 columns):"
head -5 "$OUTPUT_DIR/combined_TPM.txt" | cut -f1-5