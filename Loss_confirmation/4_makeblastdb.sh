# follow the steps to make DNA and miRNA blast database
for i in *_1.fastq.gz
do
f=${i%_1.fastq.gz}_2.fastq.gz
zcat $i $f |sed  -n '1~4s/^@/>/p;2~4p' >> SRR5055303-6.blastDB.fa
rm $i $f
done
for i in SRR5055303-6.blastDB.fa
do
makeblastdb -in $i -out $i -dbtype nucl -max_file_sz 3000000000
done
