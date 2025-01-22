### Mapping and database
for i in *.fastq.gz
do
SRA=`echo $i|sed 's/\.fastq.gz//g'`
minimap2 -t 32 -ax map-pb GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.NC_056661.1.fna $i | samtools sort -O BAM - > "$SRA".bam
done
samtools merge packbio_merged.bam *.bam
samtools sort packbio_merged.bam -o packbio_merged.sorted.bam
samtools index packbio_merged.sorted.bam

for i in *.fastq.gz
do
SRA=`echo $i|sed 's/\.fastq.gz//g'`
minimap2 -t 64 -ax map-pb GCF_003254725.2_ASM325472v2_genomic.NC_064256.1.fna $i | samtools sort -O BAM - > "$SRA".bam
done

parallel-fastq-dump --sra-id SRR25817557 --threads 16 --outdir . --split-files --gzip 


for i in *_1.fastq.gz
do
f=${i%_1.fastq.gz}_2.fastq.gz
zcat $i $f |sed  -n '1~4s/^@/>/p;2~4p' >  ${i%_1.fastq.gz}.fa
makeblastdb -in  ${i%_1.fastq.gz}.fa -out  ${i%_1.fastq.gz}.fa -dbtype nucl
done


### Canis_lupus_dingo

transcriptID=NM_017654.4
efetch -db nuccore -id $transcriptID -format fasta_cds_na|cut -f1,2 -d'_'|sed 's/>lcl|/>/g' > "$transcriptID".fa
python make_codon_position_bed.py $transcriptID
grep -w "$transcriptID" /mnt/disk4/BUDDHA/SAMD9-9L/CDK6-HEPACAM2/chain-Canis_lupus_dingo_NC_064256.1/TOGA_SAMD9/inact_mut_data.txt > "$transcriptID".inact_mut_data.txt

cat "$transcriptID".inact_mut_data.txt |cut -f3-6|grep -v "GENE:"|awk 'NF == 4 && $1!=0 {print "Exon"$1":codon"$2":"$3":"$4}' > "$transcriptID".events.bed
bash generate_codon_bed_for_cds.sh $transcriptID

for i in `cat "$transcriptID".events.bed`; do pos=`echo $i|cut -f2 -d':'`; grep -w "$pos" "$transcriptID".cds_codons.bed|sed "s/$pos/$i/g"|awk '{print $0"\t"$2"\t"$3"\t"255","0","0}'  >> "$transcriptID".codon_positions.events.bed; done



# BLASTn in DNAseqDB
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db /media/lokdeep/sdf/BUDDHA/SAMD9_loss_carnivore/Canis_lupus_dingo/Illumina_SRA/SRR25817557.fa -out "$transcriptID".blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query "$transcriptID".fa
sed -i "s/Query_1/$transcriptID/g" "$transcriptID".blastn.DNAseqDB.sam
samtools view -bhS "$transcriptID".blastn.DNAseqDB.sam > "$transcriptID".blastn.DNAseqDB.bam
samtools sort "$transcriptID".blastn.DNAseqDB.bam -o "$transcriptID".blastn.DNAseqDB.sorted.bam
samtools index "$transcriptID".blastn.DNAseqDB.sorted.bam

# create IGV report
create_report "$transcriptID".codon_positions.events.bed --standalone --fasta "$transcriptID".fa --tracks "$transcriptID".blastn.DNAseqDB.sorted.bam "$transcriptID".codon_positions.events.bed --output "$transcriptID".BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track
mkdir out_Canis_lupus_dingo
mv NM_017654.4.* out_Canis_lupus_dingo/

### Panthera_tigris
transcriptID=NM_017654.4
efetch -db nuccore -id $transcriptID -format fasta_cds_na|cut -f1,2 -d'_'|sed 's/>lcl|/>/g' > "$transcriptID".fa
python make_codon_position_bed.py $transcriptID
grep -w "$transcriptID" /mnt/disk4/BUDDHA/SAMD9-9L/CDK6-HEPACAM2/chain-Panthera_tigris_NC_056661.1/TOGA_SAMD9/inact_mut_data.txt > "$transcriptID".inact_mut_data.txt

cat "$transcriptID".inact_mut_data.txt |cut -f3-6|grep -v "GENE:"|awk 'NF == 4 && $1!=0 {print "Exon"$1":codon"$2":"$3":"$4}' > "$transcriptID".events.bed
bash generate_codon_bed_for_cds.sh $transcriptID

for i in `cat "$transcriptID".events.bed`; do pos=`echo $i|cut -f2 -d':'`; grep -w "$pos" "$transcriptID".cds_codons.bed|sed "s/$pos/$i/g"|awk '{print $0"\t"$2"\t"$3"\t"255","0","0}'  >> "$transcriptID".codon_positions.events.bed; done



# BLASTn in DNAseqDB
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db /media/lokdeep/sdf/BUDDHA/SAMD9_loss_carnivore/Panthera_tigris/Illumina_db/SRR17636806.fa -out "$transcriptID".blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query "$transcriptID".fa
sed -i "s/Query_1/$transcriptID/g" "$transcriptID".blastn.DNAseqDB.sam
samtools view -bhS "$transcriptID".blastn.DNAseqDB.sam > "$transcriptID".blastn.DNAseqDB.bam
samtools sort "$transcriptID".blastn.DNAseqDB.bam -o "$transcriptID".blastn.DNAseqDB.sorted.bam
samtools index "$transcriptID".blastn.DNAseqDB.sorted.bam

# create IGV report
create_report "$transcriptID".codon_positions.events.bed --standalone --fasta "$transcriptID".fa --tracks "$transcriptID".blastn.DNAseqDB.sorted.bam "$transcriptID".codon_positions.events.bed --output "$transcriptID".BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track




############################# LONG read verification #########################################

## Canis_lupus_dingo
# chr14:18,304,090-18,491,305
klumpy scan_alignments --alignment_map SRR11348531.bam --threads 32
klumpy find_gaps --fasta GCF_003254725.2_ASM325472v2_genomic.NC_064256.1.fna
klumpy alignment_plot --alignment_map SRR11348531.bam --reference NC_064256.1 --candidates SRR11348531_Candidate_Regions.tsv --min_len 3000 --window_size 10000 --window_step 5000 --color red --vertical_line_gaps --vertical_line_klumps --format svg --leftbound 18304090 --rightbound 18491305 --annotation GCF_003254725.2_ASM325472v2_genomic.gtf --gap_file GCF_003254725.2_ASM325472v2_genomic.NC_064256.1_gaps.tsv --width 1536 --height 864
#Parsing SRR11348531.bam
#A total of 1328 alignment records will be drawn
#Plotting..
#Parsing GCF_003254725.2_ASM325472v2_genomic.gtf
#Found the following genes: LOC112674926 LOC112674812 LOC112674782 LOC112674921 HEPACAM2 LOC125752462

## Panthera_tigris
# chrA2:95,559,116-95,687,258
klumpy scan_alignments --alignment_map packbio_merged.sorted.bam --threads 32
klumpy find_gaps --fasta GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.NC_056661.1.fna
klumpy alignment_plot --alignment_map packbio_merged.sorted.bam --reference NC_056661.1 --candidates packbio_merged.sorted_Candidate_Regions.tsv --min_len 3000 --window_size 10000 --window_step 5000 --color red --vertical_line_gaps --vertical_line_klumps --format png --leftbound 95559116 --rightbound 95687258 --annotation GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.gtf --gap_file GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.NC_056661.1_gaps.tsv --width 1536 --height 864
#Parsing packbio_merged.sorted.bam
#A total of 419 alignment records will be drawn
#Plotting..
#Parsing GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.gtf
#Found the following genes: LOC102949850 LOC102949552 SAMD9L LOC122233561 HEPACAM2
