#!/bin/bash

# Function to process the fasta file
process_fasta_file() {
    local input_file=$1
    local output_dir=$2
    local base_name=$3
    local suffix=$4
    local blast_executable="blastn"
    # Define output files
    local new_output_file="${output_dir}/${base_name}_${suffix}_results.txt"
    local temp_output_file="${new_output_file}.tmp"

    # Run BLAST for the entire file
    $blast_executable -query "$input_file" -db nt -out "$temp_output_file" -outfmt "6 qseqid sseqid qstart qend sstart send evalue bitscore pident qcovs" -max_target_seqs 3 -remote

    # Filter results to keep only the top 3 matches per query
    awk '{
        count[$1]++
        if (count[$1] <= 3) {
            print $0
        }
    }' "$temp_output_file" > "$new_output_file"

    # Remove the temporary file after processing
    rm "$temp_output_file"

    # Print a message indicating the file has been processed
    echo "Completed processing: $input_file. Result saved to $new_output_file"
}

# Main script logic
base_dir="./Q12/results" #change to input directive with barcodes (can be in multiple subdirectories)
output_base_dir="./BLAST"  # Change this to the desired output base directory
# Create the output directory if it doesn't exist
mkdir -p "$output_base_dir"

# Find all directories that might contain fasta files
for dir in $(find "$base_dir" -type d); do
    # Process files ending with _consensussequences.fasta in each directory
    for file in "$dir"/*_consensussequences.fasta; do
        if [[ -f "$file" ]]; then
            # Extract the base name without extension and suffix
            base_name=$(basename "$file" | sed -E 's/(_consensussequences)\.fasta//')
            # Count underscores in the base name
            underscore_count=$(echo "$base_name" | awk -F_ '{print NF-1}')
            # Process only files with exactly one underscore before the suffix
            if (( underscore_count == 1 )); then
                if [[ "$file" == *_consensussequences.fasta ]]; then
                    # Call the function to process the consensussequences file
                    process_fasta_file "$file" "$output_base_dir" "$base_name" "consensussequences"
                fi
            fi
        fi
    done
done
