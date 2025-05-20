#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/372/415/GCF_028372415.1_mMacEug1.pri_v2/GCF_028372415.1_mMacEug1.pri_v2_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/372/415/GCF_028372415.1_mMacEug1.pri_v2/GCF_028372415.1_mMacEug1.pri_v2_genomic.gtf.gz

#gunzip GCF_028372415.1_mMacEug1.pr*

#!/bin/bash
for i in `tail -n+2 SRA.lst`
do
parallel-fastq-dump --sra-id $i --threads 64 --outdir . --split-files --gzip
done

## for indexing
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir . --genomeFastaFiles *.fna --sjdbGTFfile *.gtf
## for mapping
for i in *_1.fastq.gz
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
STAR --runThreadN 16 --outSAMtype BAM SortedByCoordinate --genomeDir . --readFilesIn $i $j --readFilesCommand zcat --outFileNamePrefix "$k"_
samtools index "$k"_*.bam
done


rm -r chr* *.tab sjdbInfo.txt genomeParameters.txt Genome SA SAindex *.out *.gz *_STAR*
