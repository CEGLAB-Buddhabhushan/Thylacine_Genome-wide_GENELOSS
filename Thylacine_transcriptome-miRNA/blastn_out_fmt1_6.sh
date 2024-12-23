for i in Sarcophilus_harrisii* Antechinus_flavipes_*
do
blastn -task blastn -evalue 0.05 -max_target_seqs 10000 -db SRR23147611-6.miRNAblastDB.fa -out blastn."$i".SRR23147611-6.miRNAblastDB.1out -num_threads 64 -outfmt 1 -query $i
blastn -task blastn -evalue 0.05 -max_target_seqs 10000 -db SRR23147611-6.miRNAblastDB.fa -out blastn."$i".SRR23147611-6.miRNAblastDB.6out -num_threads 64 -outfmt 6 -query $i
done
