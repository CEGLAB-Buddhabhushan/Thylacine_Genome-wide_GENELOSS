for gene in CUZD1 HSD17B13  SAMD9L  VWA7
do
cd $gene
for query in *.fa
do
blastn -task blastn -evalue 0.01 -max_target_seqs 10000 -db /media/lokdeep/sdf/BUDDHA/Thylacine/out/Thylacine_genome_SRA_db/GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa -out blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.out -num_threads 64 -outfmt "7 std qseq sseq stitle" -query $query

mview -in blast -moltype dna -out fasta blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.out |sed 's/^.*SRR/>SRR/g;s/[0-9]://g'|sed 's/^.*VAHE/>VAHE/g;s/[0-9]://g'|sed 's/^.*CM/>CM/g;s/[0-9]://g'|awk -F' ' '{print $1}' > blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.fasta
mview -in fasta -html head -find 'TGA:TAA:TAG' -css on  blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.fasta > blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.html

done
cd ..
echo $gene
done
