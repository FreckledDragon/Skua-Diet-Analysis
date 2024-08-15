#!/bin/bash
# Directory containing the .fastq files
input_dir="./Q12" 
output_dir="./Q12/results"
mkdir -p "$output_dir"  # Create the output directory if it doesn't exist
# Loop through each .fastq file in the directory
for file in "$input_dir"/*.fastq; do
	#check if file exists to avoid running on non-matching patterns
	if [[ -f "$file" ]]; then
  		# Extract the base filename without extension
 		filename=$(basename "$file" .fastq)
		#Debug: print filename being processed
echo "Processing file: $filename"
 	 	# Count the number of lines that start with '@', i.e., the reads
  		line_count=$(grep -c '^@' "$file")
		#Check if count is >5 (amplicon_sorter doesn't work otherwise)
		if [[ "$line_count" -lt 5 ]]; then
			echo "Skipping $filename as line count is less than 5"
			continue
		fi
		#Get max reads
 		maxr_value=$((line_count * 100))
		#Extract file name after underscore
		new_filename="${filename#*_}"
		#Construct output file path
		output_file="$output_dir/$new_filename"
  		#Run the amplicon_sorter with the calculated maxr value
python3 amplicon_sorter.py -i "$file" -ra -maxr "$maxr_value" -np 10 -o "$output_file" -min 100
	fi
done
echo "Processing complete. "
