#!/usr/bin/env bash

# rank_mhc_samples.sh - Rank samples by MHC expression for organoid selection

INPUT_FILE="${1:-mhc_expression_all_samples.tsv}"
OUTPUT_DIR="${OUTPUT_DIR:-mhc_ranking_results}"
TOP_N="${TOP_N:-10}"

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
info(){ echo "[$(ts)] [INFO ] $*"; }
warn(){ echo "[$(ts)] [WARN ] $*" >&2; }
die(){ echo "[$(ts)] [ERROR] $*" >&2; exit 1; }

# Check input file exists
[[ -f "$INPUT_FILE" ]] || die "Input file not found: $INPUT_FILE"

# Create output directory
mkdir -p "$OUTPUT_DIR"

info "Ranking MHC expression from: $INPUT_FILE"
info "Output directory: $OUTPUT_DIR"
info "Selecting top $TOP_N samples"

# Calculate total MHC expression per sample
info "Calculating total MHC expression per sample..."
awk -F'\t' '
    NR == 1 { next }  # Skip header
    {
        sample = $1
        gene_name = $2
        tpm = ($4 != "" && $4 != "NA") ? $4 : 0
        expected_count = ($5 != "" && $5 != "NA") ? $5 : 0
        
        # Accumulate TPM and expected counts
        sample_tpm[sample] += tpm
        sample_expected[sample] += expected_count
        gene_count[sample]++
        
        # Track individual gene expressions for detailed output
        sample_genes[sample, gene_name] = tpm "\t" expected_count "\t" $6
    }
    END {
        # Output ranking by total TPM
        print "sample\ttotal_TPM\ttotal_expected_count\tnum_mhc_genes\tavg_TPM_per_gene"
        for (sample in sample_tpm) {
            total_tpm = sample_tpm[sample]
            total_expected = sample_expected[sample]
            num_genes = gene_count[sample]
            avg_tpm = (num_genes > 0) ? total_tpm / num_genes : 0
            print sample "\t" total_tpm "\t" total_expected "\t" num_genes "\t" avg_tpm
        }
    }
' "$INPUT_FILE" | sort -k2,2nr > "$OUTPUT_DIR/sample_ranking_by_total_TPM.tsv"

# Also rank by average TPM per gene
info "Calculating average MHC expression per sample..."
awk -F'\t' '
    NR == 1 { next }
    {
        sample = $1
        tpm = ($4 != "" && $4 != "NA") ? $4 : 0
        sample_tpm[sample] += tpm
        gene_count[sample]++
    }
    END {
        print "sample\tavg_TPM_per_gene\ttotal_TPM\tnum_mhc_genes"
        for (sample in sample_tpm) {
            total_tpm = sample_tpm[sample]
            num_genes = gene_count[sample]
            avg_tpm = (num_genes > 0) ? total_tpm / num_genes : 0
            print sample "\t" avg_tpm "\t" total_tpm "\t" num_genes
        }
    }
' "$INPUT_FILE" | sort -k2,2nr > "$OUTPUT_DIR/sample_ranking_by_avg_TPM.tsv"

# Create detailed profiles for top samples
info "Creating detailed profiles for top $TOP_N samples..."

# Get top samples by total TPM
top_samples=($(tail -n +2 "$OUTPUT_DIR/sample_ranking_by_total_TPM.tsv" | head -n "$TOP_N" | cut -f1))

info "Top $TOP_N samples by total MHC expression: ${top_samples[*]}"

# Create detailed expression profiles
{
    echo -e "rank\tsample\tgene_name\tTPM\texpected_count\tSTAR_counts"
    rank=1
    for sample in "${top_samples[@]}"; do
        awk -F'\t' -v sample="$sample" -v rank="$rank" '
            $1 == sample {
                print rank "\t" $1 "\t" $2 "\t" $4 "\t" $5 "\t" $6
            }
        ' "$INPUT_FILE" | sort -k4,4nr  # Sort by TPM within each sample
        ((rank++))
    done
} > "$OUTPUT_DIR/top_${TOP_N}_detailed_profiles.tsv"

# Create summary statistics
info "Generating summary statistics..."
{
    echo "=== MHC EXPRESSION RANKING SUMMARY ==="
    echo "Generated: $(date)"
    echo "Input file: $INPUT_FILE"
    echo "Total samples analyzed: $(tail -n +2 "$INPUT_FILE" | cut -f1 | sort -u | wc -l)"
    echo "Total MHC genes found: $(tail -n +2 "$INPUT_FILE" | cut -f2 | sort -u | wc -l)"
    echo ""
    echo "TOP $TOP_N SAMPLES FOR ORGANOID CULTURE (by total MHC TPM):"
    echo "Rank | Sample | Total_TPM | Avg_TPM | Genes | Recommendation"
    echo "-----|--------|-----------|---------|-------|---------------"
    
    rank=1
    while IFS=$'\t' read -r sample total_tpm total_expected num_genes avg_tpm; do
        if (( rank <= TOP_N )); then
            # Add recommendation based on expression levels
            if (( $(echo "$total_tpm > 200" | bc -l) )); then
                rec="Excellent - High expression"
            elif (( $(echo "$total_tpm > 100" | bc -l) )); then
                rec="Good - Moderate expression"  
            elif (( $(echo "$total_tpm > 50" | bc -l) )); then
                rec="Fair - Low expression"
            else
                rec="Poor - Very low expression"
            fi
            
            printf "%4d | %8s | %9.2f | %7.2f | %5d | %s\n" "$rank" "$sample" "$total_tpm" "$avg_tpm" "$num_genes" "$rec"
        fi
        ((rank++))
    done < <(tail -n +2 "$OUTPUT_DIR/sample_ranking_by_total_TPM.tsv")
    
    echo ""
    echo "KEY MHC GENES TO MONITOR:"
    echo "- HLA-A, HLA-B, HLA-C: Classical MHC Class I (most important)"
    echo "- HLA-DRB1, HLA-DQB1: Key MHC Class II genes"
    echo "- MICA, MICB: NK cell ligands"
    
    echo ""
    echo "FILES GENERATED:"
    echo "- sample_ranking_by_total_TPM.tsv: Ranked by total MHC expression"
    echo "- sample_ranking_by_avg_TPM.tsv: Ranked by average expression per gene"
    echo "- top_${TOP_N}_detailed_profiles.tsv: Detailed gene-by-gene profiles"
    echo "- summary_report.txt: This summary report"
    
} > "$OUTPUT_DIR/summary_report.txt"

# Create a simple top 10 list for quick reference
{
    echo "TOP $TOP_N SAMPLES FOR ORGANOID CULTURE"
    echo "Sample_ID"
    tail -n +2 "$OUTPUT_DIR/sample_ranking_by_total_TPM.tsv" | head -n "$TOP_N" | cut -f1
} > "$OUTPUT_DIR/top_${TOP_N}_samples.txt"

# Create gene-specific rankings (which samples have highest expression of key genes)
info "Creating gene-specific rankings..."
{
    echo -e "gene_name\ttop_sample\tmax_TPM\ttop_3_samples"
    for gene in HLA-A HLA-B HLA-C HLA-DRB1 MICA MICB; do
        top_info=$(awk -F'\t' -v gene="$gene" '
            $2 == gene && $4 != "NA" { 
                samples[NR] = $1 "\t" $4
                tpm[NR] = $4
            }
            END {
                # Find max TPM
                max_tpm = 0
                max_sample = ""
                for (i in tpm) {
                    if (tpm[i] > max_tpm) {
                        max_tpm = tpm[i]
                        split(samples[i], parts, "\t")
                        max_sample = parts[1]
                    }
                }
                
                # Get top 3 samples for this gene
                cmd = "echo"
                for (i in samples) {
                    cmd = cmd " \"" samples[i] "\""
                }
                cmd = cmd " | tr \" \" \"\\n\" | sort -k2,2nr | head -3 | cut -f1 | tr \"\\n\" \",\""
                cmd | getline top3
                close(cmd)
                gsub(/,$/, "", top3)  # Remove trailing comma
                
                print gene "\t" max_sample "\t" max_tpm "\t" top3
            }
        ' "$INPUT_FILE")
        echo "$top_info"
    done
} > "$OUTPUT_DIR/gene_specific_rankings.tsv"

info "=== RANKING COMPLETE ==="
info "Output files created in: $OUTPUT_DIR/"
info "Key files:"
info "  - top_${TOP_N}_samples.txt: Simple list of top samples"
info "  - sample_ranking_by_total_TPM.tsv: Full ranking with scores"
info "  - top_${TOP_N}_detailed_profiles.tsv: Detailed expression profiles"
info "  - summary_report.txt: Human-readable summary"
info "  - gene_specific_rankings.tsv: Best samples for each key gene"

# Display quick results
echo ""
echo "QUICK RESULTS - Top $TOP_N samples for organoid culture:"
cat "$OUTPUT_DIR/top_${TOP_N}_samples.txt"

echo ""
echo "Summary statistics:"
head -15 "$OUTPUT_DIR/summary_report.txt" | tail -10