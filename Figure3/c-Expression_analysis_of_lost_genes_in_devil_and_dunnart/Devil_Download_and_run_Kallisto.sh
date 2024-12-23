grep "Illumina HiSeq" SRA.table.tsv | awk -F'\t' '{split($10, a, "."); print a[1]"-"$7}' > SRA.dw.lst
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

## download cds fasta 
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/635/505/GCF_902635505.1_mSarHar1.11/GCF_902635505.1_mSarHar1.11_cds_from_genomic.fna.gz
#gunzip GCF_902635505.1_mSarHar1.11_cds_from_genomic.fna.gz
##change the header of CDS file
#sed 's/^.*gene=/>/' GCF_902635505.1_mSarHar1.11_cds_from_genomic.fna | cut -f1,4 -d']'|sed 's/\[protein_id=//g;s/\]//g;s/ /-/g' > Reformated.GCF_902635505.1_mSarHar1.11_cds_from_genomic.fna
#for i in `cat protein.lst`;do id=`esearch -db protein -query $i|  efetch -format gp|grep "DBSOURCE"|sed 's/DBSOURCE    REFSEQ: accession //g'`; sed -i "s/$i/$i-$id/g" Reformated.GCF_902635505.1_mSarHar1.11_cds_from_genomic.fna; echo $i-$id;done

kallisto index --make-unique -i transcripts.idx Reformated.GCF_902635505.1_mSarHar1.11_cds_from_genomic.fna
for i in `ls -1 *_1.fastq.gz`
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
kallisto quant -t 32 -i transcripts.idx -o "$k".kallisto_out -b 100 $i $j
done
