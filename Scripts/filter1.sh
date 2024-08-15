#!/bin/bash
# Define input and output directories
input_dir="./barcodes"
output_dir="./Q12"
# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
# Compress all .fastq files in the input directory
for file in "$input_dir"/*.fastq; do
    if [ -e "$file" ]; then
        gzip "$file"
fi
done
# Process each .fastq.gz file using NanoFilt and save the output uncompressed
echo "Processing .fastq.gz files in $input_dir"
for file in "$input_dir"/*.fastq.gz; do
    if [ -e "$file" ]; then
        og_base_name=$(basename "$file" .fastq.gz)
        # Extract part of the name after the first underscore
        base_name="${og_base_name#*_}"
        gunzip -c "$file" | NanoFilt -q 12 -l 100 --maxlength 700 > "$output_dir/$base_name.fastq"
        echo "Processed $file"
    else
        echo "No .fastq.gz files found in $input_dir"
    fi
done
echo "All files have been processed."
