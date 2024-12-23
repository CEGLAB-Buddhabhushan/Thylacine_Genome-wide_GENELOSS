query="Sarcophilus_harrisii_VWA7_XM_003768930.3_CDS.fa"
cut -f2,3,5 outfmt6/blastn."$query".*|awk '{print ">"$1"\n"$3}'|sed 's/-//g' > All.fa
makeblastdb -in All.fa -out All.fa -dbtype nucl
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db All.fa -out blastn."$query".All.fa.mview.out -num_threads 64 -outfmt "7 std qseq sseq stitle" -query $query

mview -in blast -moltype dna -out fasta  blastn."$query".All.fa.mview.out |sed 's/^.*SRR/>SRR/g;s/[0-9]://g'|sed 's/^.*VAHE/>VAHE/g;s/[0-9]://g'|sed 's/^.*CM/>CM/g;s/[0-9]://g'|awk -F' ' '{print $1}' > blastn."$query".All.fa.mview.out.fasta
mview -in fasta -html head -find 'TGA:TAA:TAG' -css on  blastn."$query".All.fa.mview.out.fasta > blastn."$query".All.fa.mview.out.fasta.html
