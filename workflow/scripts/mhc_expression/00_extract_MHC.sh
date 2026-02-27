#!/usr/bin/env bash

# Simple working batch script - uses direct gene ID matching

STAR_COL="${STAR_COL:-2}"
EXCLUDE="${EXCLUDE:-log logs tmp temp}"
COHORT_LONG="mhc_expression_all_samples.tsv"
SKIP_REPORT="mhc_batch_skip_report.tsv"

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
info(){ echo "[$(ts)] [INFO ] $*"; }
warn(){ echo "[$(ts)] [WARN ] $*" >&2; }

# Initialize output files
echo -e "sample\tgene_name\tgene_id\tTPM\texpected_count\tSTAR_counts" > "$COHORT_LONG"
echo -e "sample\treason\tfound_rsem\tfound_star\tfound_gtf" > "$SKIP_REPORT"

process_one() {
    local d="$1"
    local sample="$(basename "$d")"
    
    info "Processing sample: $sample"
    
    # Check for required files
    local rsem_file="$d/${sample}_rsem.genes.results"
    local star_file="$d/${sample}_ReadsPerGene.out.tab"
    
    local missing=()
    [[ -f "$rsem_file" ]] || missing+=("rsem")
    [[ -f "$star_file" ]] || missing+=("star")
    
    if (( ${#missing[@]} > 0 )); then
        echo -e "${sample}\tMissing: ${missing[*]}\t${rsem_file}\t${star_file}\tNA" >> "$SKIP_REPORT"
        warn "  Skipping $sample - missing: ${missing[*]}"
        return 1
    fi
    
    # Extract MHC genes using known gene IDs - direct approach
    awk -F'\t' -v sample="$sample" -v star_col="$STAR_COL" '
        BEGIN {
            # Define MHC gene mappings
            mhc_map["ENSG00000206503"] = "HLA-A"
            mhc_map["ENSG00000234745"] = "HLA-B"
            mhc_map["ENSG00000204525"] = "HLA-C"
            mhc_map["ENSG00000204520"] = "MICA"
            mhc_map["ENSG00000204516"] = "MICB"
            mhc_map["ENSG00000204257"] = "HLA-DMA"
            mhc_map["ENSG00000242574"] = "HLA-DMB"
            mhc_map["ENSG00000204252"] = "HLA-DOA"
            mhc_map["ENSG00000241106"] = "HLA-DOB"
            mhc_map["ENSG00000231389"] = "HLA-DPA1"
            mhc_map["ENSG00000231461"] = "HLA-DPA2"
            mhc_map["ENSG00000237398"] = "HLA-DPA3"
            mhc_map["ENSG00000223865"] = "HLA-DPB1"
            mhc_map["ENSG00000224557"] = "HLA-DPB2"
            mhc_map["ENSG00000196735"] = "HLA-DQA1"
            mhc_map["ENSG00000237541"] = "HLA-DQA2"
            mhc_map["ENSG00000179344"] = "HLA-DQB1"
            mhc_map["ENSG00000232629"] = "HLA-DQB2"
            mhc_map["ENSG00000226030"] = "HLA-DQB3"
            mhc_map["ENSG00000204287"] = "HLA-DRA"
            mhc_map["ENSG00000196126"] = "HLA-DRB1"
            mhc_map["ENSG00000198502"] = "HLA-DRB5"
            mhc_map["ENSG00000229391"] = "HLA-DRB6"
            mhc_map["ENSG00000196301"] = "HLA-DRB9"
            mhc_map["ENSG00000204592"] = "HLA-E"
            mhc_map["ENSG00000204642"] = "HLA-F"
            mhc_map["ENSG00000204632"] = "HLA-G"
            mhc_map["ENSG00000206341"] = "HLA-H"
            mhc_map["ENSG00000204622"] = "HLA-J"
            mhc_map["ENSG00000230795"] = "HLA-K"
            mhc_map["ENSG00000243753"] = "HLA-L"
            mhc_map["ENSG00000224372"] = "HLA-N"
            mhc_map["ENSG00000261548"] = "HLA-P"
            mhc_map["ENSG00000225851"] = "HLA-S"
            mhc_map["ENSG00000231130"] = "HLA-T"
            mhc_map["ENSG00000228078"] = "HLA-U"
            mhc_map["ENSG00000181126"] = "HLA-V"
            mhc_map["ENSG00000235290"] = "HLA-W"
            mhc_map["ENSG00000235301"] = "HLA-Z"
        }
        
        # Process RSEM file
        FILENAME ~ /rsem/ && FNR > 1 {
            gid = $1
            gid_base = gid
            gsub(/\.[0-9]+$/, "", gid_base)
            
            if (gid_base in mhc_map) {
                gene_name = mhc_map[gid_base]
                tpm = $6
                expected_count = $5
                rsem_data[gid_base] = sample "\t" gene_name "\t" gid "\t" tpm "\t" expected_count
            }
        }
        
        # Process STAR file  
        FILENAME ~ /ReadsPerGene/ && $1 !~ /^N_/ {
            gid = $1
            gid_base = gid
            gsub(/\.[0-9]+$/, "", gid_base)
            
            if (gid_base in rsem_data) {
                star_count = (star_col <= NF ? $star_col : "NA")
                print rsem_data[gid_base] "\t" star_count
            }
        }
    ' "$rsem_file" "$star_file" >> "$COHORT_LONG"
    
    info "  Completed $sample"
    return 0
}

# Main processing loop
info "Starting batch processing in: $(pwd)"
processed=0
skipped=0

# Get list of directories  
dirs=()
for dir in */; do
    dir_name="${dir%/}"
    [[ -d "$dir_name" ]] || continue
    
    # Check exclusions
    excluded=false
    for x in $EXCLUDE; do
        if [[ "$dir_name" == "$x" ]]; then
            excluded=true
            break
        fi
    done
    [[ "$excluded" == "true" ]] && continue
    
    dirs+=("$dir_name")
done

info "Found ${#dirs[@]} directories to process"

# Process each directory
for d in "${dirs[@]}"; do
    if process_one "$d"; then
        ((processed++))
    else
        ((skipped++))
    fi
    
    if (( processed % 10 == 0 )); then
        info "Progress checkpoint: $processed processed, $skipped skipped"
    fi
done

info "=== BATCH COMPLETE ==="
info "Successfully processed: $processed samples"
info "Skipped: $skipped samples"
info "Total data lines: $(($(wc -l < "$COHORT_LONG") - 1))"

# Show some results
if (( processed > 0 )); then
    info "Sample results (first 10 lines):"
    head -11 "$COHORT_LONG"
    
    # Build matrices
    info "Building matrices..."
    
    # TPM matrix
    awk -F'\t' 'BEGIN{OFS="\t"} 
        NR==1{next} 
        {
            sample=$1; gene=$2; tpm=$4
            samples[sample]=1; genes[gene]=1; data[gene,sample]=tpm
        } 
        END{
            # Print header
            printf "gene_name"
            for(s in samples) printf OFS s
            print ""
            
            # Print data rows
            for(g in genes) {
                printf g
                for(s in samples) printf OFS (data[g,s] ? data[g,s] : "NA")
                print ""
            }
        }' "$COHORT_LONG" > "mhc_TPM_matrix.tsv"
    
    info "Created mhc_TPM_matrix.tsv"
    info "Matrix dimensions: $(head -1 mhc_TPM_matrix.tsv | awk -F'\t' '{print NF-1}') samples x $(tail -n +2 mhc_TPM_matrix.tsv | wc -l) genes"
fi

info "All done!"