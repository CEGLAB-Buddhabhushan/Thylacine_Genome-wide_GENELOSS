#!/bin/bash
transcriptID=$1
# Input FASTA file
fasta_file="$transcriptID".fa

# Extract sequence name and length using samtools faidx
sequence_name=$(grep "^>" $fasta_file | sed 's/>//')
sequence_length=$(faidx $fasta_file -i chromsizes | awk '{print $2}')

# Define the output BED6 file
output_file="${sequence_name}.cds_codons.bed"

# Initialize start position
start=0

# Open the output file
echo -n "" > $output_file

# Generate the BED6 file
for (( i=1; i<=$sequence_length; i+=3 ))
do
    end=$((start + 3))
    codon_name="codon$(( (i + 3) / 3 ))"
    echo -e "${sequence_name}\t${start}\t${end}\t${codon_name}\t0\t+" >> $output_file
    start=$((start + 3))
done

echo "BED6 file generated: $output_file"
