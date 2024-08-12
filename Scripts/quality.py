import os
import csv
from Bio import SeqIO

# Define the input directory containing FASTQ files and the output CSV file
input_dir = "./fastq_files/"  # Replace with your directory containing FASTQ files
output_file = "./fastq_files/scores.csv"  # Output CSV file name

# Create a CSV file to store the results
with open(output_file, mode="w", newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["File", "Sequence ID", "Average Quality Score", "Sequence Length"])

    # Loop over all FASTQ files in the input directory
    for fastq_file in os.listdir(input_dir):
        if fastq_file.endswith(".fastq"):
            input_path = os.path.join(input_dir, fastq_file)
            
            # Parse the FASTQ file and calculate average quality scores and sequence lengths
            for record in SeqIO.parse(input_path, "fastq"):
                quality_scores = record.letter_annotations["phred_quality"]
                average_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
                sequence_length = len(record.seq)
                csv_writer.writerow([fastq_file, record.id, average_quality, sequence_length])

print(f"Quality scores and sequence lengths have been written to {output_file}")
