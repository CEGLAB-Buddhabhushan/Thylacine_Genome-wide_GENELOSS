cut -f2 HSD17B11/outfmt6/blastn.Sarcophilus_harrisii_HSD17B11_CDS.fa.SRR5055303-6.blastDB.6out|sort -u > diff_btw_HSD17B11_HSD17B13/HSD17B11_CDS.read.list
cut -f2 outfmt6/blastn.Sarcophilus_harrisii_HSD17B13_Exons.fa.SRR5055303-6.blastDB.6out|sort -u > diff_btw_HSD17B11_HSD17B13/HSD17B13_Exons.read.list

cd diff_btw_HSD17B11_HSD17B13/
comm -23 HSD17B13_Exons.read.list HSD17B11_CDS.read.list|sort -u > HSD17B13.uniq_read.list

for i in SRR5055303 SRR5055304 SRR5055305 SRR5055306
do
grep "$i" HSD17B13.uniq_read.list > "$i".lst
seqtk subseq out/"$i"_1.fastq.gz "$i".lst > "$i"_1.HSD17B13.fq
seqtk subseq out/"$i"_2.fastq.gz "$i".lst > "$i"_2.HSD17B13.fq
done

sed -n '1~4s/^@/>/p;2~4p' *.HSD17B13.fq > SRR5055303-6.HSD17B13.fa
makeblastdb -in SRR5055303-6.HSD17B13.fa -out SRR5055303-6.HSD17B13.fa -dbtype nucl

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.HSD17B13.fa -out blastn.Sarcophilus_harrisii_HSD17B13_Exons.fa.SRR5055303-6.HSD17B13.fa.1out -num_threads 64 -outfmt 1 -query ../Sarcophilus_harrisii_HSD17B13_Exons.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.HSD17B13.fa -out blastn.Sarcophilus_harrisii_HSD17B13_Exons.fa.SRR5055303-6.HSD17B13.fa.4out -num_threads 64 -outfmt 4 -query ../Sarcophilus_harrisii_HSD17B13_Exons.fa

blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.HSD17B13.fa -out blastn.Sarcophilus_harrisii_HSD17B13_Exons_with_Intron.fa.SRR5055303-6.HSD17B13.fa.1out -num_threads 64 -outfmt 1 -query ../Sarcophilus_harrisii_HSD17B13_Exons_with_Intron.fa
blastn -task blastn -evalue 0.05 -max_target_seqs 5000 -db SRR5055303-6.HSD17B13.fa -out blastn.Sarcophilus_harrisii_HSD17B13_Exons_with_Intron.fa.SRR5055303-6.HSD17B13.fa.4out -num_threads 64 -outfmt 4 -query ../Sarcophilus_harrisii_HSD17B13_Exons_with_Intron.fa

mkdir igv_report
cp igv_report

blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db ../SRR5055303-6.HSD17B13.fa -out XM_003774508.4.blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query XM_003774508.4.fa
sed -i "s/Query_1/XM_003774508.4/g" XM_003774508.4.blastn.DNAseqDB.sam
samtools view -bhS XM_003774508.4.blastn.DNAseqDB.sam > XM_003774508.4.blastn.DNAseqDB.bam
samtools sort XM_003774508.4.blastn.DNAseqDB.bam -o XM_003774508.4.blastn.DNAseqDB.sorted.bam
samtools index XM_003774508.4.blastn.DNAseqDB.sorted.bam


# create IGV report
create_report XM_003774508.4.codon_positions.events.bed --standalone --fasta XM_003774508.4.fa --tracks XM_003774508.4.blastn.DNAseqDB.sorted.bam XM_003774508.4.codon_positions.events.bed --output XM_003774508.4.BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track


blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db ../SRR5055303-6.HSD17B13.fa -out Exon1_with_intron.blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query Exon1_with_intron.fa
sed -i "s/Query_1/Exon1/g" Exon1_with_intron.blastn.DNAseqDB.sam
samtools view -bhS Exon1_with_intron.blastn.DNAseqDB.sam > Exon1_with_intron.blastn.DNAseqDB.bam
samtools sort Exon1_with_intron.blastn.DNAseqDB.bam -o Exon1_with_intron.blastn.DNAseqDB.sorted.bam
samtools index Exon1_with_intron.blastn.DNAseqDB.sorted.bam


# create IGV report
create_report Exon1_with_intron.codon_positions.events.bed --standalone --fasta Exon1_with_intron.fa --tracks Exon1_with_intron.blastn.DNAseqDB.sorted.bam Exon1_with_intron.codon_positions.events.bed --output Exon1_with_intron.BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track

