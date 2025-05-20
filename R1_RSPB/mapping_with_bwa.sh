bwa index GCF_902635505.1_mSarHar1.11_genomic.fna
for i in *_1.fastq.gz; do j=`echo $i|sed 's/_1/_2/g'`; k=`echo $i|sed 's/_1.fastq.gz//g'`
bwa mem -t 32 GCF_902635505.1_mSarHar1.11_genomic.fna $i $j |samtools sort -@ 32 -O BAM -o "$k".bwa.bam
samtools index -c "$k".bwa.bam
done
samtools merge -@ 16 SRR5055303-6.merged.bam *.bwa.bam


for i in *_1.fastq.gz
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
samtools index -c "$k".bwa.bam
samtools flagstat -@ 16 "$k".bwa.bam > "$k".bwa.flagstat.out
done
for i in *_1.fastq.gz
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
samtools index -@ 16 -c "$k".bwa.bam
samtools flagstat -@ 16 "$k".bwa.bam > "$k".bwa.flagstat.out
done
samtools merge -@ 16 SRR5055303-6.merged.bam *.bwa.bam
rm SRR5055303-6.merged.bam
samtools merge -@ 32 SRR5055303-6.bwa.merged.bam *.bwa.bam
samtools sort SRR5055303-6.bwa.merged.bam -@ 24 -o SRR5055303-6.bwa.sorted.bam
samtools index -c -@ 16 SRR5055303-6.bwa.sorted.bam
rm SRR5055303-6.bwa.merged.bam 





for i in *_1.fastq.gz
do
j=`echo $i|sed 's/_1/_2/g'`
k=`echo $i|sed 's/_1.fastq.gz//g'`
minimap2 -t 32 -ax sr GCF_902635505.1_mSarHar1.11_genomic.fna $i $j | samtools sort -@ 32 -O BAM -o "$k".minimap2.bam
samtools index -@ 16 -c "$k".minimap2.bam
samtools flagstat -@ 16 "$k".minimap2.bam > "$k".minimap2.flagstat.out
done
samtools merge -@ 16 SRR5055303-6.minimap2.merged.bam *.minimap2.bam
samtools sort SRR5055303-6.minimap2.merged.bam -@ 16 -o SRR5055303-6.minimap2.sorted.bam
samtools index -@ 16 -c SRR5055303-6.minimap2.sorted.bam
rm SRR5055303-6.minimap2.merged.bam *.minimap2.bam



###stapy
/home/buddha/software/stampy/stampy.py -G SarHar GCF_902635505.1_mSarHar1.11_genomic.fna
/home/buddha/software/stampy/stampy.py -g SarHar -H SarHar
/home/buddha/software/stampy/stampy.py -g SarHar -h SarHar -t 36 --bamkeepgoodreads -M SRR5055303-6.bwa.sorted.bam |samtools sort -@ 16 -O BAM -o SRR5055303-6.bwa.sorted.bamkeepgoodreads_stampy.bam

##########Monodelphis_domestica 
# blastn databse
zcat *.hifi_reads.fastq.gz|sed -n '1~4s/^@/>/p;2~4p' > Monodelphis_domestica.pacbio.BLASTdb.fa && makeblastdb -in Monodelphis_domestica.pacbio.BLASTdb.fa -out Monodelphis_domestica.pacbio.BLASTdb.fa -dbtype nucl  -max_file_sz 3G
### pacbio_hifi assembly verification
bwa index GCF_027887165.1_mMonDom1.pri_genomic.fna
for i in *.hifi_reads.fastq.gz
do
k=`echo $i|cut -f1 -d'.'`
bwa mem -t 64 -x pacbio GCF_027887165.1_mMonDom1.pri_genomic.fna $i  |samtools sort -@ 64 -O BAM -o "$k".bwa.pacbio.bam
samtools index -c "$k".bwa.pacbio.bam
done
samtools merge -@ 32 MonDom1.merged.bam *.bwa.pacbio.bam
samtools sort MonDom1.merged.bam  -@ 48 -o MonDom1.bwa.pacbio.sorted.bam
samtools index -c -@ 48 MonDom1.bwa.pacbio.sorted.bam
rm MonDom1.merged.bam *.bwa.pacbio.bam *.hifi_reads.fastq.gz


