#!/bin/bash
input_file="./BLAST/result.txt"
output_file="./BLAST/final.txt"

#sort sequences by unique barcode contig and unique species by highest identical score
awk '
{
    key = $1
    if (!(key in max) || $9 > max[key]) {
        max[key] = $9
        lines[key] = $0
    }
}
END {
    for (key in lines) {
        print lines[key]
    }
}
' "$input_file" > "$output_file"


