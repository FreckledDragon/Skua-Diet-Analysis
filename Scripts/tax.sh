#!/bin/bash
# Define paths to files
blast_output="./BLAST/filtered_blast.txt"
seq_ids_file="./sequence_ids.txt"
taxonomy_info_file="./taxonomy_info.txt"
combined_output="./BLAST/result.txt"
# Create or empty the output files if they already exist
> "$seq_ids_file"
> "$taxonomy_info_file"
> "$combined_output"
# Extract sseqid from BLAST results without removing duplicates
awk -F'\t' '{split($2, a, "|"); if (length(a) >= 4) print a[4]}' "$blast_output" > "$seq_ids_file"
# Extract protein
for protein in $(cat "$seq_ids_file"); do
	result=$(esearch -db nucleotide -query "$protein" | efetch -format docsum | xtract -pattern DocumentSummary -element Organism)
	if [ -z "$result" ]; then
		echo "NA"
	else
		echo "$result"
	fi
done > "$taxonomy_info_file"
paste "$blast_output" "$seq_ids_file" "$taxonomy_info_file" > "$combined_output"

