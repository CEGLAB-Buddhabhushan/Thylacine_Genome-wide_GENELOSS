STAR --runThreadN 16 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Myrmecobius_fasciatus.VWA7.JAJPUD010000745.1.fa --genomeSAindexNbases 8
## for mapping
for i in *_1.fastq.gz
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
STAR --runThreadN 16 --outSAMtype BAM SortedByCoordinate --genomeDir . --readFilesIn $i $j --readFilesCommand zcat --outFileNamePrefix "$k"_ --limitBAMsortRAM 1031735638
samtools index "$k"_*.bam
done
