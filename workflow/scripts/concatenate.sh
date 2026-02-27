#!/bin/bash
#SBATCH -t 36:00:00
#SBATCH --mem=15G
#SBATCH --account=hansengroup
#SBATCH -J concatenate
#SBATCH -p all
#SBATCH -c 4
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS

set -euo pipefail

# Usage: sbatch concatenate.sh input_manifest.csv output_manifest.tsv
# Input manifest format: sample_name,group_name,mate1,mate2
#   - group_name non-empty: samples with the same group_name get concatenated
#   - group_name empty:     sample stays separate (falls back to sample_name)
INPUT_MANIFEST="${1:-260225_manifest_replicates.csv}"
OUTPUT_MANIFEST="${2:-manifest_final.tsv}"
OUTPUT_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process/fastq"
COMBINED_MANIFEST="manifest_combined.csv"

mkdir -p "$OUTPUT_DIR"

echo -e "sample\tread1\tread2\tsource_sample_names" > "$OUTPUT_MANIFEST"
echo "sample_name,mate1,mate2,group" > "$COMBINED_MANIFEST"

# Build lookup: group_name -> source_group (folder under .../Data/RNAseq/)
declare -A group_to_source
while IFS=',' read -r sample_name group_name mate1 mate2; do
    [[ "$sample_name" == "sample_name" ]] && continue
    grp="${group_name:-$sample_name}"
    source_group=$(echo "$mate1" | sed -n 's|.*/cluster/projects/livingbank/Data/RNAseq/\([^/]*\)/.*|\1|p')
    if [[ -n "$source_group" && -n "$grp" ]]; then
        [[ -z "${group_to_source[$grp]:-}" ]] && group_to_source[$grp]="$source_group"
    fi
done < <(tr -d '\r' < "$INPUT_MANIFEST")

# Build a sorted temp file: group_name<TAB>mate1<TAB>mate2<TAB>sample_name
# sort -k1,1 groups by group_name while preserving within-group CSV order (stable sort)
TMP_FILE="$(mktemp)"
tr -d '\r' < "$INPUT_MANIFEST" | awk -F ',' 'NR==1 { next } {
    sample_name = $1
    group_name  = $2
    mate1       = $3
    mate2       = $4
    grp = (group_name != "") ? group_name : sample_name
    print grp "\t" mate1 "\t" mate2 "\t" sample_name
}' | sort -k1,1 > "$TMP_FILE"

current_group=""
declare -a r1_files=()
declare -a r2_files=()
declare -a source_names=()

flush_group() {
    [[ -z "$current_group" || ${#r1_files[@]} -eq 0 ]] && return

    echo "Processing: $current_group (${#r1_files[@]} file(s))" >&2

    local first_r1="${r1_files[0]}"
    local ext=""
    case "$first_r1" in
        *.fastq.gz) ext=".fastq.gz" ;;
        *.fq.gz)    ext=".fq.gz"    ;;
        *.fastq)    ext=".fastq"    ;;
        *.fq)       ext=".fq"       ;;
    esac

    local combined_r1="${OUTPUT_DIR}/${current_group}_R1${ext}"
    local combined_r2="${OUTPUT_DIR}/${current_group}_R2${ext}"

    : > "$combined_r1"
    for f in "${r1_files[@]}"; do
        if [[ -r "$f" ]]; then
            cat "$f" >> "$combined_r1"
        else
            echo "WARNING: R1 not readable: '$f'" >&2
        fi
    done

    : > "$combined_r2"
    for f in "${r2_files[@]}"; do
        if [[ -r "$f" ]]; then
            cat "$f" >> "$combined_r2"
        else
            echo "WARNING: R2 not readable: '$f'" >&2
        fi
    done

    if [[ -s "$combined_r1" && -s "$combined_r2" ]]; then
        local source_join
        source_join="$(IFS=','; echo "${source_names[*]}")"
        echo -e "${current_group}\t${combined_r1}\t${combined_r2}\t${source_join}" >> "$OUTPUT_MANIFEST"
        local source_group="${group_to_source[$current_group]:-unknown}"
        echo "${current_group},${combined_r1},${combined_r2},${source_group}" >> "$COMBINED_MANIFEST"
        echo "Done: $current_group" >&2
    else
        echo "WARNING: no readable files for '$current_group', skipping" >&2
        rm -f "$combined_r1" "$combined_r2"
    fi

    r1_files=()
    r2_files=()
    source_names=()
}

while IFS=$'\t' read -r grp mate1 mate2 sample_name; do
    if [[ "$grp" != "$current_group" ]]; then
        flush_group
        current_group="$grp"
    fi

    if [[ ! -r "$mate1" || ! -r "$mate2" ]]; then
        echo "WARNING: skipping '$sample_name' — file(s) unreadable (mate1='$mate1', mate2='$mate2')" >&2
        continue
    fi

    r1_files+=("$mate1")
    r2_files+=("$mate2")
    source_names+=("$sample_name")
done < "$TMP_FILE"

flush_group

rm -f "$TMP_FILE"

echo "" >&2
echo "Output manifest:   $OUTPUT_MANIFEST" >&2
echo "Combined manifest: $COMBINED_MANIFEST" >&2
