4#!/bin/bash
# Define directory with files
input_dir="./BLAST"
# Define output file
output_file="./filtered_blast.txt"
# Clear the output file to avoid appending previous data
> "$output_file"
# Create output directory
mkdir -p "$input_dir/res"
# Debugging: Print paths to check values
echo "Input directory: $input_dir"
echo "Output file: $output_file"
# Loop through all .txt files in the specified directory
for file in "$input_dir"/*.txt; do
  if [ -f "$file" ]; then
    echo "Processing file: $file"
    # Apply the awk command and append the results to the output file
    awk -F'\t' '$9 > 97 {print}' "$file" >> "$output_file"
  fi
done
echo "Processing complete. Results saved in $output_file."
