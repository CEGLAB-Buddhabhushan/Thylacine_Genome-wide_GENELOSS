#To validate inactivating mutations, we extracted the genomic sequence around the mutation with a 25-bp flank and used megablast (parameters match score −1, mismatch scores −2, gap costs linear, expectation value threshold 10) to search unassembled sequencing reads stored in the TRACE and Sequence Read Archives (supplementary table 2, Supplementary Material online). 

#sensitive blastn runs (word size = 7) to search read data of


blastn -task megablast \
-query mutation_flank.fa \
-db read_db \
-out megablast_results.txt \
-evalue 10 \
-reward -1 \
-penalty -2 \
-gapopen 0 \
-gapextend 2 \
-outfmt 6


blastn -task blastn \
-query mutation_flank.fa \
-db read_db \
-word_size 7 \
-evalue 10 \
-out sensitive_blastn_results.txt \
-outfmt 6


#https://pmc.ncbi.nlm.nih.gov/articles/PMC5716171/ BlastN in “megablast” mode with default parameters
blastn -task megablast \
-query mutation_flank.fa \
-db read_db \
-out megablast_default_results.txt


###################################### BLAST: Rechecking of Gene Loss Events ####################################################
for gene in CUZD1 HSD17B13 SAMD9L VWA7
do
cd $gene/BLAST
for transcriptID in `ls -1`
do
mkdir "$gene"-Rechecking_with_blast
cd "$gene"-Rechecking_with_blast
mkdir $transcriptID
cd $transcriptID
efetch -db nuccore -id $transcriptID -format fasta_cds_na|cut -f1,2 -d'_'|sed 's/>lcl|/>/g' > "$transcriptID".fa
blastn -task blastn -evalue 0.05  -max_target_seqs 5000 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out "$transcriptID".blastn.evalue0.05.sam -num_threads 64 -outfmt '17 SQ'  -query "$transcriptID".fa
blastn -task blastn -evalue 0.01  -max_target_seqs 5000 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out "$transcriptID".blastn.evalue0.01.sam -num_threads 64 -outfmt '17 SQ'  -query "$transcriptID".fa
blastn -task blastn -evalue 10 -word_size 7  -max_target_seqs 5000 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out "$transcriptID".blastn.evalue10.word_size7.sam -num_threads 64 -outfmt '17 SQ'  -query "$transcriptID".fa
blastn -task megablast -evalue 10 -reward 1 -penalty -2 -max_target_seqs 5000 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out "$transcriptID".megablast.evalue10.match1.mismatch-2.sam -num_threads 64 -outfmt '17 SQ'  -query "$transcriptID".fa
blastn -task megablast -evalue 10 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out "$transcriptID".megablast.evalue10.sam -num_threads 64 -outfmt '17 SQ'  -query "$transcriptID".fa
blastn -task dc-megablast -evalue 10 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out "$transcriptID".dc-megablast.evalue10.sam -num_threads 64 -outfmt '17 SQ'  -query "$transcriptID".fa

#convert sam to bam 
for sam in *.sam
do
sed -i "s/Query_1/$transcriptID/g" $sam
samtools view -bhS $sam > ${sam%.sam}.bam
samtools sort ${sam%.sam}.bam -o ${sam%.sam}.sorted.bam
samtools index ${sam%.sam}.sorted.bam
done
bams=`ls -1 *.sorted.bam|tr '\n' ' '`
cp ../../$transcriptID/"$transcriptID".codon_positions.events.bed .
##create IGV-report 
create_report "$transcriptID".codon_positions.events.bed --standalone --fasta "$transcriptID".fa --tracks $bams "$transcriptID".codon_positions.events.bed --output "$transcriptID".BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track
cd ../..
done
cd /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/LOST_GENES/
echo $gene
done


##################################### Mapping: Rechecking of Gene Loss Events ################################################# 
for gene in CUZD1 HSD17B13 SAMD9L VWA7
do
cd $gene/
mkdir "$gene"-Rechecking_with_mapping
cd "$gene"-Rechecking_with_mapping
for transcriptID in `ls -1 ../Genome_mapping/|cut -f1,2 -d'.'|cut -f1,2 -d'_'|sort -u`
do
##create IGV-report 
cp ../Genome_mapping/"$transcriptID"_codon_positions.events.bed .

create_report "$transcriptID"_codon_positions.events.bed --standalone --fasta /home/buddha/work/R1_Thylacine/Thylacine_SRA_DNA/GCF_902635505.1_mSarHar1.11_genomic.fna --tracks /home/buddha/work/R1_Thylacine/Thylacine_SRA_DNA/SRR5055303-6.bwa.sorted.bam /home/buddha/work/R1_Thylacine/Thylacine_SRA_DNA/SRR5055303-6.minimap2.sorted.bam /home/buddha/work/R1_Thylacine/Thylacine_SRA_DNA/SRR5055303-6.bowtie2.sorted.bam "$transcriptID"_codon_positions.events.bed /home/buddha/work/R1_Thylacine/Thylacine_SRA_DNA/GCF_902635505.1_mSarHar1.11_genomic.gtf --output "$transcriptID".Genome_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track
done
cd ../../
echo $gene
done

