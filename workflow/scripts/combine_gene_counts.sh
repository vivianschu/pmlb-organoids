#!/bin/bash

# =============================================================================
# Extract gene detection counts from RSEM results
# =============================================================================

BASE_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
OUTPUT_DIR="${BASE_DIR}/combined_results"

mkdir -p "$OUTPUT_DIR"

OUTPUT_FILE="$OUTPUT_DIR/gene_detection_stats.tsv"

echo -e "sample\ttotal_genes\tgenes_detected_TPM0\tgenes_detected_TPM0.1\tgenes_detected_TPM1\tgenes_detected_TPM10" > "$OUTPUT_FILE"

# Find all RSEM gene results
# Files are in ${BASE_DIR}/{sample_name}/{sample_name}_rsem.genes.results
RSEM_FILES=($(find "$BASE_DIR" -maxdepth 2 -name "*_rsem.genes.results" -type f | sort))

for file in "${RSEM_FILES[@]}"; do
    sample=$(basename "$file" "_rsem.genes.results")
    echo "Processing: $sample"

    # Compute TPM from expected_count (col 5) and effective_length (col 4):
    #   RPK = expected_count / (effective_length / 1000)
    #   TPM = (RPK / sum(RPK)) * 1e6
    # Then count genes at different TPM thresholds in a single pass.
    read -r total_genes detected_0 detected_0_1 detected_1 detected_10 < <(
        awk -F'\t' '
            NR==1 { next }
            {
                count[NR] = $5
                eff_len[NR] = $4
                rpk[NR] = (eff_len[NR] > 0) ? count[NR] / (eff_len[NR] / 1000) : 0
                sum_rpk += rpk[NR]
                total++
            }
            END {
                d0=0; d01=0; d1=0; d10=0
                for (i=2; i<=NR; i++) {
                    tpm = (sum_rpk > 0) ? (rpk[i] / sum_rpk) * 1e6 : 0
                    if (tpm >  0)   d0++
                    if (tpm >= 0.1) d01++
                    if (tpm >= 1)   d1++
                    if (tpm >= 10)  d10++
                }
                print total, d0, d01, d1, d10
            }
        ' "$file"
    )

    echo -e "${sample}\t${total_genes}\t${detected_0}\t${detected_0_1}\t${detected_1}\t${detected_10}" >> "$OUTPUT_FILE"
done

echo ""
echo "Done! Output: $OUTPUT_FILE"
head "$OUTPUT_FILE" | column -t