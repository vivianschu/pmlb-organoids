#!/usr/bin/env python3
"""
Parse the RNA-seq multisample CSV file and output a simple TSV format
that bash can easily read.

Output format (TSV):
sample_name    mate1_path    mate2_path    group_name

For samples with multiple replicates, outputs separate lines with _A, _B, _C suffixes.
The group_name column equals the sample_name (grouping is determined by the
manifest CSV's group_name column, handled downstream by concatenate.sh).

This parser handles the non-standard CSV format where bracketed Python lists
contain commas but are not properly quoted.
"""

import sys
import ast


def parse_line(line):
    """
    Parse a single line from the CSV file.
    Format: organoid_id,num_replicates,[path1, path2, ...],[path1, path2, ...]
    
    The challenge is that the lists contain commas but aren't quoted in the CSV.
    """
    line = line.strip()
    if not line:
        return None
    
    # Strategy: Find the first '[' which marks the start of R1 paths
    # Everything before that (minus trailing comma) is "organoid_id,num_replicates"
    
    first_bracket = line.find('[')
    if first_bracket == -1:
        return None
    
    # Get the prefix part (organoid_id,num_replicates,)
    prefix = line[:first_bracket].rstrip(',')
    parts = prefix.split(',')
    
    if len(parts) < 2:
        return None
    
    organoid_id = parts[0]
    try:
        num_replicates = int(parts[1])
    except ValueError:
        return None
    
    # Now extract the two bracketed lists
    # Find the pattern: [...],[...]
    rest = line[first_bracket:]
    
    # Find matching brackets for R1
    bracket_count = 0
    r1_end = -1
    for i, char in enumerate(rest):
        if char == '[':
            bracket_count += 1
        elif char == ']':
            bracket_count -= 1
            if bracket_count == 0:
                r1_end = i
                break
    
    if r1_end == -1:
        return None
    
    r1_str = rest[:r1_end + 1]
    
    # R2 starts after "],["
    remaining = rest[r1_end + 1:].lstrip(',')
    r2_str = remaining
    
    # Parse the lists using ast.literal_eval
    try:
        r1_paths = ast.literal_eval(r1_str)
        r2_paths = ast.literal_eval(r2_str)
    except (ValueError, SyntaxError) as e:
        print(f"Warning: Could not parse paths for {organoid_id}: {e}", file=sys.stderr)
        return None
    
    # Ensure they're lists
    if not isinstance(r1_paths, list):
        r1_paths = [r1_paths] if r1_paths else []
    if not isinstance(r2_paths, list):
        r2_paths = [r2_paths] if r2_paths else []
    
    return {
        'organoid_id': organoid_id,
        'num_replicates': num_replicates,
        'r1_paths': r1_paths,
        'r2_paths': r2_paths
    }



def parse_csv(csv_file, output_file=None):
    """Parse the CSV and output sample information."""

    samples = []

    with open(csv_file, 'r') as f:
        # Skip header
        f.readline()

        for line_num, line in enumerate(f, start=2):
            parsed = parse_line(line)

            if parsed is None:
                continue

            organoid_id = parsed['organoid_id']
            num_replicates = parsed['num_replicates']
            r1_paths = parsed['r1_paths']
            r2_paths = parsed['r2_paths']

            # Skip samples with 0 replicates
            if num_replicates == 0:
                continue

            # Validate that we have the expected number of paths
            if len(r1_paths) != num_replicates or len(r2_paths) != num_replicates:
                print(f"Warning: Path count mismatch for {organoid_id} (line {line_num}): "
                      f"expected {num_replicates}, got R1={len(r1_paths)}, R2={len(r2_paths)}",
                      file=sys.stderr)

            # Create entries for each replicate
            for i in range(min(len(r1_paths), len(r2_paths))):
                # Add suffix for multiple replicates (_A, _B, _C, etc.)
                if num_replicates > 1:
                    sample_name = f"{organoid_id}_{chr(65 + i)}"  # A, B, C, ...
                else:
                    sample_name = organoid_id

                samples.append({
                    'name': sample_name,
                    'mate1': r1_paths[i],
                    'mate2': r2_paths[i],
                    'group_name': sample_name
                })
    
    # Output
    output = sys.stdout if output_file is None else open(output_file, 'w')

    # Write header
    output.write("sample_name\tmate1\tmate2\tgroup_name\n")

    # Write samples
    for sample in samples:
        output.write(f"{sample['name']}\t{sample['mate1']}\t{sample['mate2']}\t{sample['group_name']}\n")
    
    if output_file:
        output.close()
    
    # Print summary to stderr
    print(f"Total samples (including replicates): {len(samples)}", file=sys.stderr)
    
    return samples


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input_csv> [output_tsv]", file=sys.stderr)
        sys.exit(1)

    csv_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    parse_csv(csv_file, output_file)


if __name__ == '__main__':
    main()