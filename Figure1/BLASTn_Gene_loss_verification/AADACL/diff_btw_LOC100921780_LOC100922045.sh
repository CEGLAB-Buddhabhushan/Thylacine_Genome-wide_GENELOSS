mkdir diff_btw_LOC100921780_LOC100922045
cut -f2 outfmt6/blastn.Sarcophilus_harrisii_LOC100921780_ENSSHAT00000036342.1.fa.SRR5055303-6.blastDB.6out|sort -u > diff_btw_LOC100921780_LOC100922045/LOC100921780_Exons.read.list
cut -f2 outfmt6/blastn.Sarcophilus_harrisii_LOC100922045_XM_003766079.2_CDS.fa.SRR5055303-6.blastDB.6out|sort -u > diff_btw_LOC100921780_LOC100922045/LOC100922045_Exons.read.list

cd diff_btw_LOC100921780_LOC100922045
comm -23 LOC100922045_Exons.read.list LOC100921780_Exons.read.list|sort -u > LOC100922045.uniq_read.list

for i in SRR5055303 SRR5055304 SRR5055305 SRR5055306
do
grep "$i" LOC100922045.uniq_read.list > "$i".lst
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA_output/BLASTn_Gene_loss_verification/HSD17B13/diff_btw_HSD17B11_HSD17B13/out/"$i"_1.fastq.gz "$i".lst > "$i"_1.LOC100922045.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA_output/BLASTn_Gene_loss_verification/HSD17B13/diff_btw_HSD17B11_HSD17B13/out/"$i"_2.fastq.gz "$i".lst > "$i"_2.LOC100922045.fq
done

sed -n '1~4s/^@/>/p;2~4p' *.LOC100922045.fq > SRR5055303-6.LOC100922045.fa
makeblastdb -in SRR5055303-6.LOC100922045.fa -out SRR5055303-6.LOC100922045.fa -dbtype nucl

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100922045.fa -out blastn.Sarcophilus_harrisii_LOC100922045_XM_003766079.2_CDS.fa.SRR5055303-6.LOC100922045.fa.1out -num_threads 64 -outfmt 1 -query ../Sarcophilus_harrisii_LOC100922045_XM_003766079.2_CDS.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100922045.fa -out blastn.Sarcophilus_harrisii_LOC100922045_XM_003766079.2_CDS.fa.SRR5055303-6.LOC100922045.fa.4out -num_threads 64 -outfmt 4 -query ../Sarcophilus_harrisii_LOC100922045_XM_003766079.2_CDS.fa

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100922045.fa -out blastn.Sarcophilus_harrisii_LOC100922045_ENSSHAT00000020538.2_Exons.fa.SRR5055303-6.LOC100922045.fa.1out -num_threads 64 -outfmt 1 -query ../Sarcophilus_harrisii_LOC100922045_ENSSHAT00000020538.2_Exons.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100922045.fa -out blastn.Sarcophilus_harrisii_LOC100922045_ENSSHAT00000020538.2_Exons.fa.SRR5055303-6.LOC100922045.fa.4out -num_threads 64 -outfmt 4 -query ../Sarcophilus_harrisii_LOC100922045_ENSSHAT00000020538.2_Exons.fa

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100922045.fa -out blastn.Antechinus_flavipes_LOC127556989_XM_051990490.1_CDS.fa.SRR5055303-6.LOC100922045.fa.1out -num_threads 64 -outfmt 1 -query ../Antechinus_flavipes_LOC127556989_XM_051990490.1_CDS.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100922045.fa -out blastn.Antechinus_flavipes_LOC127556989_XM_051990490.1_CDS.fa.SRR5055303-6.LOC100922045.fa.4out -num_threads 64 -outfmt 4 -query ../Antechinus_flavipes_LOC127556989_XM_051990490.1_CDS.fa


mkdir igv_report
cd igv_report

blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db ../SRR5055303-6.LOC100922045.fa -out XM_003766079.2.blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query XM_003766079.2.fa
sed -i "s/Query_1/XM_003766079.2/g" XM_003766079.2.blastn.DNAseqDB.sam
samtools view -bhS XM_003766079.2.blastn.DNAseqDB.sam > XM_003766079.2.blastn.DNAseqDB.bam
samtools sort XM_003766079.2.blastn.DNAseqDB.bam -o XM_003766079.2.blastn.DNAseqDB.sorted.bam
samtools index XM_003766079.2.blastn.DNAseqDB.sorted.bam


# create IGV report
create_report XM_003766079.2.codon_positions.events.bed --standalone --fasta XM_003766079.2.fa --tracks XM_003766079.2.blastn.DNAseqDB.sorted.bam XM_003766079.2.codon_positions.events.bed --output XM_003766079.2.BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track

