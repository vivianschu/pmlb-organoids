#!/bin/bash

# =============================================================================
# Combine STAR Log.final.out files from multiple samples
# =============================================================================

# Directory containing all sample folders
BASE_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"  # Change this to your actual path
OUTPUT_DIR="${BASE_DIR}/combined_results"

mkdir -p "$OUTPUT_DIR"

OUTPUT_FILE="$OUTPUT_DIR/STAR_mapping_stats.tsv"

# Find all Log.final.out files
# Files are in ${BASE_DIR}/{sample_name}/{sample_name}_Log.final.out
LOG_FILES=($(find "$BASE_DIR" -maxdepth 2 -name "*_Log.final.out" -type f | sort))

if [ ${#LOG_FILES[@]} -eq 0 ]; then
    echo "No STAR Log.final.out files found!"
    exit 1
fi

echo "Found ${#LOG_FILES[@]} samples"

# Write header
cat > "$OUTPUT_FILE" << 'EOF'
sample	total_reads	avg_input_read_length	uniquely_mapped_reads	uniquely_mapped_percent	avg_mapped_length	num_splices	num_splices_annotated	num_splices_GT_AG	num_splices_GC_AG	num_splices_AT_AC	num_splices_non_canonical	mismatch_rate	deletion_rate	deletion_avg_length	insertion_rate	insertion_avg_length	multimapped_reads	multimapped_percent	multimapped_toomany_reads	multimapped_toomany_percent	unmapped_mismatches_percent	unmapped_tooshort_percent	unmapped_other_percent	chimeric_reads	chimeric_percent
EOF

# Function to extract value from log file
extract_value() {
    local file="$1"
    local pattern="$2"
    grep "$pattern" "$file" | awk -F'|' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' | tr -d '%'
}

# Process each log file
for file in "${LOG_FILES[@]}"; do
    # Extract sample name from file path
    sample=$(basename "$file" "_Log.final.out")
    
    echo "Processing: $sample"
    
    # Extract all metrics
    total_reads=$(extract_value "$file" "Number of input reads")
    avg_input_length=$(extract_value "$file" "Average input read length")
    unique_reads=$(extract_value "$file" "Uniquely mapped reads number")
    unique_percent=$(extract_value "$file" "Uniquely mapped reads %")
    avg_mapped_length=$(extract_value "$file" "Average mapped length")
    num_splices=$(extract_value "$file" "Number of splices: Total")
    num_splices_annotated=$(extract_value "$file" "Number of splices: Annotated")
    num_splices_GT_AG=$(extract_value "$file" "Number of splices: GT/AG")
    num_splices_GC_AG=$(extract_value "$file" "Number of splices: GC/AG")
    num_splices_AT_AC=$(extract_value "$file" "Number of splices: AT/AC")
    num_splices_noncanon=$(extract_value "$file" "Number of splices: Non-canonical")
    mismatch_rate=$(extract_value "$file" "Mismatch rate per base")
    deletion_rate=$(extract_value "$file" "Deletion rate per base")
    deletion_avg_len=$(extract_value "$file" "Deletion average length")
    insertion_rate=$(extract_value "$file" "Insertion rate per base")
    insertion_avg_len=$(extract_value "$file" "Insertion average length")
    multi_reads=$(extract_value "$file" "Number of reads mapped to multiple loci")
    multi_percent=$(extract_value "$file" "% of reads mapped to multiple loci")
    multi_toomany_reads=$(extract_value "$file" "Number of reads mapped to too many loci")
    multi_toomany_percent=$(extract_value "$file" "% of reads mapped to too many loci")
    unmapped_mismatch=$(extract_value "$file" "% of reads unmapped: too many mismatches")
    unmapped_short=$(extract_value "$file" "% of reads unmapped: too short")
    unmapped_other=$(extract_value "$file" "% of reads unmapped: other")
    chimeric_reads=$(extract_value "$file" "Number of chimeric reads")
    chimeric_percent=$(extract_value "$file" "% of chimeric reads")
    
    # Write to output
    echo -e "${sample}\t${total_reads}\t${avg_input_length}\t${unique_reads}\t${unique_percent}\t${avg_mapped_length}\t${num_splices}\t${num_splices_annotated}\t${num_splices_GT_AG}\t${num_splices_GC_AG}\t${num_splices_AT_AC}\t${num_splices_noncanon}\t${mismatch_rate}\t${deletion_rate}\t${deletion_avg_len}\t${insertion_rate}\t${insertion_avg_len}\t${multi_reads}\t${multi_percent}\t${multi_toomany_reads}\t${multi_toomany_percent}\t${unmapped_mismatch}\t${unmapped_short}\t${unmapped_other}\t${chimeric_reads}\t${chimeric_percent}" >> "$OUTPUT_FILE"
done

echo ""
echo "Done! Output file: $OUTPUT_FILE"
echo ""

# Create a summary stats file (simplified version for quick QC)
SUMMARY_FILE="$OUTPUT_DIR/STAR_mapping_summary.tsv"

echo -e "sample\ttotal_reads\tuniquely_mapped_percent\tmultimapped_percent\tunmapped_percent" > "$SUMMARY_FILE"

tail -n +2 "$OUTPUT_FILE" | while IFS=$'\t' read -r sample total_reads avg_input unique_reads unique_pct rest; do
    # Extract multi and unmapped percentages
    multi_pct=$(echo "$rest" | cut -f14)
    unmapped_short=$(echo "$rest" | cut -f18)
    unmapped_other=$(echo "$rest" | cut -f19)
    
    # Calculate total unmapped (approximately)
    unmapped_pct=$(awk "BEGIN {printf \"%.2f\", 100 - $unique_pct - $multi_pct}")
    
    echo -e "${sample}\t${total_reads}\t${unique_pct}\t${multi_pct}\t${unmapped_pct}"
done >> "$SUMMARY_FILE"

echo "Summary file: $SUMMARY_FILE"
echo ""

# Preview
echo "Preview of mapping stats (first 5 samples):"
head -6 "$OUTPUT_FILE" | cut -f1-5 | column -t

echo ""
echo "Key columns for Figure 1:"
echo "  - total_reads: for sequencing depth"
echo "  - uniquely_mapped_percent: for mapping rate"
echo "  - Use with gene counts for detection rate"