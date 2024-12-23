mkdir diff_btw_LOC100924681_LOC100924419
cut -f2 outfmt6/blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.SRR5055303-6.blastDB.6out|sort -u > diff_btw_LOC100924681_LOC100924419/LOC100924681_CDS.read.list
cut -f2 outfmt6/blastn.Sarcophilus_harrisii_LOC100924419_XM_003771945.4_CDS.fa.SRR5055303-6.blastDB.6out|sort -u > diff_btw_LOC100924681_LOC100924419/LOC100924419_CDS.read.list

cd diff_btw_LOC100924681_LOC100924419
comm -23 LOC100924681_CDS.read.list LOC100924419_CDS.read.list|sort -u > LOC100924681.uniq_read.list

for i in SRR5055303 SRR5055304 SRR5055305 SRR5055306
do
grep "$i" LOC100924681.uniq_read.list > "$i".lst
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA_output/BLASTn_Gene_loss_verification/HSD17B13/diff_btw_HSD17B11_HSD17B13/out/"$i"_1.fastq.gz "$i".lst > "$i"_1.LOC100924681.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA_output/BLASTn_Gene_loss_verification/HSD17B13/diff_btw_HSD17B11_HSD17B13/out/"$i"_2.fastq.gz "$i".lst > "$i"_2.LOC100924681.fq
done

sed -n '1~4s/^@/>/p;2~4p' *.LOC100924681.fq > SRR5055303-6.LOC100924681.fa
makeblastdb -in SRR5055303-6.LOC100924681.fa -out SRR5055303-6.LOC100924681.fa -dbtype nucl

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100924681.fa -out blastn.Antechinus_flavipes_SAMD9L_XM_052001347.1_CDS.fa.SRR5055303-6.LOC100924681.fa.1out -num_threads 64 -outfmt 1 -query ../Antechinus_flavipes_SAMD9L_XM_052001347.1_CDS.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100924681.fa -out blastn.Antechinus_flavipes_SAMD9L_XM_052001347.1_CDS.fa.SRR5055303-6.LOC100924681.fa.4out -num_threads 64 -outfmt 4 -query ../Antechinus_flavipes_SAMD9L_XM_052001347.1_CDS.fa

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100924681.fa -out blastn.Sarcophilus_harrisii_LOC100924681_XM_003771946.4_CDS.fa.SRR5055303-6.LOC100924681.fa.1out -num_threads 64 -outfmt 1 -query ../Sarcophilus_harrisii_LOC100924681_XM_003771946.4_CDS.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100924681.fa -out blastn.Sarcophilus_harrisii_LOC100924681_XM_003771946.4_CDS.fa.SRR5055303-6.LOC100924681.fa.4out -num_threads 64 -outfmt 4 -query ../Sarcophilus_harrisii_LOC100924681_XM_003771946.4_CDS.fa

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100924681.fa -out blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.SRR5055303-6.LOC100924681.fa.1out -num_threads 64 -outfmt 1 -query ../Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.LOC100924681.fa -out blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.SRR5055303-6.LOC100924681.fa.4out -num_threads 64 -outfmt 4 -query ../Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa


mkdir igv_report
cp igv_report

blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db ../SRR5055303-6.LOC100924681.fa -out XM_003771946.4.blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query XM_003771946.4.fa
sed -i "s/Query_1/XM_003771946.4/g" XM_003771946.4.blastn.DNAseqDB.sam
samtools view -bhS XM_003771946.4.blastn.DNAseqDB.sam > XM_003771946.4.blastn.DNAseqDB.bam
samtools sort XM_003771946.4.blastn.DNAseqDB.bam -o XM_003771946.4.blastn.DNAseqDB.sorted.bam
samtools index XM_003771946.4.blastn.DNAseqDB.sorted.bam


# create IGV report
create_report XM_003771946.4.codon_positions.events.bed --standalone --fasta XM_003771946.4.fa --tracks XM_003771946.4.blastn.DNAseqDB.sorted.bam XM_003771946.4.codon_positions.events.bed --output XM_003771946.4.BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track


blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db ../SRR5055303-6.LOC100924681.fa -out XM_031940019.1.blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query XM_031940019.1.fa
sed -i "s/Query_1/XM_031940019.1/g" XM_031940019.1.blastn.DNAseqDB.sam
samtools view -bhS XM_031940019.1.blastn.DNAseqDB.sam > XM_031940019.1.blastn.DNAseqDB.bam
samtools sort XM_031940019.1.blastn.DNAseqDB.bam -o XM_031940019.1.blastn.DNAseqDB.sorted.bam
samtools index XM_031940019.1.blastn.DNAseqDB.sorted.bam


# create IGV report
create_report XM_031940019.1.codon_positions.events.bed --standalone --fasta XM_031940019.1.fa --tracks XM_031940019.1.blastn.DNAseqDB.sorted.bam XM_031940019.1.codon_positions.events.bed --output XM_031940019.1.BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track


