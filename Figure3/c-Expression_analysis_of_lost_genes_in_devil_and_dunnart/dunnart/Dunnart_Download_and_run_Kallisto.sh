grep "Illumina" SRA.table.tsv |cut -f5,6|cut -f3 -d':'|sed 's/ //g'|sed 's/\t/-/g'|sort -u > SRA.dw.lst
for i in `cat SRA.dw.lst`
do
name=`echo $i|cut -f1 -d'-'`
fq1=`echo $i|cut -f2 -d'-'|cut -f1 -d';'`
fq2=`echo $i|cut -f2 -d'-'|cut -f2 -d';'`
wget $fq1
wget $fq2
fqn=`echo $i|cut -f2 -d'-'|cut -f1 -d';'|cut -f6 -d'/'`
mv "$fqn"_1.fastq.gz "$name"_"$fqn"_1.fastq.gz
mv "$fqn"_2.fastq.gz "$name"_"$fqn"_2.fastq.gz
done

## download cds fasta https://figshare.unimelb.edu.au/projects/Fat_tailed_dunnart_transcriptome_reference_files/183307
kallisto index --make-unique -i transcripts.idx DUN.CDS.fasta
for i in `ls -1 *_1.fastq.gz`
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
kallisto quant -t 32 -i transcripts.idx -o "$k".kallisto_out -b 100 $i $j
done
