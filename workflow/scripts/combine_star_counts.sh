#!/bin/bash
# usage: bash combine_star_counts.sh
#
# Combines STAR ReadsPerGene.out.tab files across all samples into a single
# integer count matrix, selecting the correct strand column per sample based
# on the strandedness inferred during alignment.
#
# Output: combined_results/combined_star_counts.txt
#   - Rows: genes (ENSG IDs), first 4 STAR summary rows excluded
#   - Columns: samples
#   - Values: integer read counts

BASE_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
OUTPUT_DIR="${BASE_DIR}/combined_results"

mkdir -p "$OUTPUT_DIR"

# Find all ReadsPerGene files
RPG_FILES=($(find "$BASE_DIR" -maxdepth 2 -name "*_ReadsPerGene.out.tab" -type f | sort))

if [ ${#RPG_FILES[@]} -eq 0 ]; then
    echo "No ReadsPerGene.out.tab files found!"
    exit 1
fi

echo "Found ${#RPG_FILES[@]} samples"

# Extract gene IDs from the first file (skip first 4 STAR summary rows)
# Col 1 is gene_id; rows 1-4 are N_unmapped, N_multimapping, N_noFeature, N_ambiguous
echo "gene_id" > "$OUTPUT_DIR/gene_ids.tmp"
tail -n +5 "${RPG_FILES[0]}" | cut -f1 >> "$OUTPUT_DIR/gene_ids.tmp"

MISSING_STRAND=0

for file in "${RPG_FILES[@]}"; do
    sample=$(basename "$file" "_ReadsPerGene.out.tab")
    sample_dir=$(dirname "$file")
    strand_file="${sample_dir}/${sample}_strandedness.txt"

    # Determine which column to use
    # STAR ReadsPerGene columns:
    #   col 2: unstranded
    #   col 3: forward strand (HTSeq -s yes)
    #   col 4: reverse strand (HTSeq -s reverse)
    if [[ -s "$strand_file" ]]; then
        strand=$(cat "$strand_file" | tr -d '[:space:]')
    else
        echo "WARNING: strandedness file not found for ${sample}, defaulting to unstranded (col 2)" >&2
        strand="no"
        MISSING_STRAND=1
    fi

    case "$strand" in
        yes)     col=3 ;;
        reverse) col=4 ;;
        no|*)    col=2 ;;
    esac

    echo "  ${sample}: strandedness='${strand}', using col ${col}"

    # Write header (sample name) + counts (skip first 4 summary rows)
    echo "$sample" > "$OUTPUT_DIR/${sample}_counts.tmp"
    tail -n +5 "$file" | cut -f"$col" >> "$OUTPUT_DIR/${sample}_counts.tmp"
done

# Paste gene_ids + all sample columns into matrix
paste "$OUTPUT_DIR/gene_ids.tmp" "$OUTPUT_DIR"/*_counts.tmp > "$OUTPUT_DIR/combined_star_counts.txt"

# Cleanup
rm -f "$OUTPUT_DIR"/*.tmp

echo ""
echo "Done!"
echo "  Output: $OUTPUT_DIR/combined_star_counts.txt"
echo "  Samples: ${#RPG_FILES[@]}"
if [[ $MISSING_STRAND -eq 1 ]]; then
    echo "  WARNING: Some samples were missing strandedness files — check output above"
fi

echo ""
echo "Preview (first 5 rows, first 5 columns):"
head -5 "$OUTPUT_DIR/combined_star_counts.txt" | cut -f1-5 | column -t