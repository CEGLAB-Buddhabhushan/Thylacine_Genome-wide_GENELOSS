#!/bin/bash
# Create header for the output TSV file
echo -e "Gene_name\tTranscriptID\tChromosome\tTranscript_length\tGC Content\tGC_Stretch\tMissing_exon\tDeleted_exon\tTOGA_status" > S.harrisii.info.tsv

# Iterate over filtered transcript IDs
for i in $(grep -vE "unassigned_transcript|XR|NR" loss_summ_data.reformated.tsv | cut -f1); do
    TranscriptID="$i"

    # Get the gene name
    Gene_name=$(awk -v var1="$TranscriptID" '$3 == "transcript" && $0 ~ "transcript_id \"" var1 "\"" { match($0, /gene_id "([^"]+)"/, arr); print arr[1] }' GCF_902635505.1_mSarHar1.11_genomic.gtf)

    # Get chromosome information
    Chromosome=$(grep "$TranscriptID" GCF_902635505.1_mSarHar1.11_genomic.gtf | awk '$3 == "transcript"' | cut -f1)

    # Extract the transcript sequence
    grep -A1 ">ref:$TranscriptID" nucleotide.reformated.fasta | sed 's/:/_/g' > "$TranscriptID.fa"

    # Get the length and GC content
    length=$(faidx "$TranscriptID.fa" -i chromsizes | cut -f2)
    GC=$(seqkit fx2tab --name --gc "$TranscriptID.fa" | cut -f2)

    # Calculate GC stretch
    GC_Stretch=$(perl /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/GC_Stretch_finder.pl "$TranscriptID.fa" | awk '{print $5}')

    # Get inactive mutation data
    grep "$TranscriptID" inact_mut_data.txt > "$TranscriptID.inact_mut_data.txt"

    # Check for missing and deleted exons
    exon_missing="No"
    if grep -q "Missing exon" "$TranscriptID.inact_mut_data.txt"; then
        exon_missing="Yes"
    fi

    exon_deleted="No"
    if grep -q "Deleted exon" "$TranscriptID.inact_mut_data.txt"; then
        exon_deleted="Yes"
    fi

    # Get TOGA status
    toga_status=$(grep -w "$TranscriptID" loss_summ_data.reformated.tsv | cut -f2)

    # Append results to the output file
    echo -e "$Gene_name\t$TranscriptID\t$Chromosome\t$length\t$GC\t$GC_Stretch\t$exon_missing\t$exon_deleted\t$toga_status" >> S.harrisii.info.tsv

    # Print TranscriptID for tracking
    echo "$TranscriptID"

    # Clean up temporary files
    rm "$TranscriptID.fa" "$TranscriptID.inact_mut_data.txt" "$TranscriptID.fa.fai" 
done

#awk 'BEGIN{FS="\t"} NF != 9 {print "Row", NR, "has", NF, "columns."}' S.harrisii.info.tsv
## Gel longest isoform
head -n1 S.harrisii.info.tsv > S.harrisii.longest_isoform.info.tsv
tail -n +2 S.harrisii.info.tsv | \
cut -f1 | \
sort -u | \
while read GENE; do
    LONGEST=$(grep -w "^$GENE" S.harrisii.info.tsv | \
              sort -t$'\t' -k4,4nr | \
              head -n1)
    printf "$LONGEST\n" >> longest.list
done
cat longest.list | while read LONGEST_ISOFOMR; do
    grep -F "$LONGEST_ISOFOMR" S.harrisii.info.tsv >> tmp.tsv
done
sort -k1,1 tmp.tsv >> S.harrisii.longest_isoform.info.tsv
rm longest.list tmp.tsv
