#!/bin/bash

# Define input and output directories
input_dir="./path"
output_dir="./path"

# Create the output directory if it doesn't exist
mkdir -p "output_dir"

# Find all .pod5 files in subdirectories of the input directory
find "$input_dir" -type f -name "*.pod5" | while read -r file; do

  # Get the base name without the path
  basename=$(basename "$file")
  # Copy the file to the output directory with the name
  cp "$file" "$output_dir/$basename"
done
