query=$1
mkdir outfmt1 outfmt4 outfmt6 outfmt7 mview
SRA_db="/media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa"
Refassem="/media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA_output/BLASTn_Gene_loss_verification/Genomes_database/GCA_007646695.1_UniMelb_Thylacine_Refassem_1_genomic.fna"
hybrid_assembly="/media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA_output/BLASTn_Gene_loss_verification/Genomes_database/GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly_genomic.fna"

#OUT fmt 6
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $SRA_db -out outfmt6/blastn."$query".SRR5055303-6.blastDB.6out -num_threads 64 -outfmt "6 qseqid sseqid qstart qend sseq" -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $Refassem -out outfmt6/blastn."$query".GCA_007646695.1_UniMelb_Thylacine_Refassem.6out -num_threads 64 -outfmt "6 qseqid sseqid qstart qend sseq" -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $hybrid_assembly -out outfmt6/blastn."$query".GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly.6out -num_threads 64 -outfmt "6 qseqid sseqid qstart qend sseq" -query $query

#OUT fmt 1
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $SRA_db -out outfmt1/blastn."$query".SRR5055303-6.blastDB.1out -num_threads 64 -outfmt 1 -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $Refassem -out outfmt1/blastn."$query".GCA_007646695.1_UniMelb_Thylacine_Refassem.1out -num_threads 64 -outfmt 1 -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $hybrid_assembly -out outfmt1/blastn."$query".GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly.1out -num_threads 64 -outfmt 1 -query $query

#OUT fmt 4
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $SRA_db -out outfmt4/blastn."$query".SRR5055303-6.blastDB.4out -num_threads 64 -outfmt 4 -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $Refassem -out outfmt4/blastn."$query".GCA_007646695.1_UniMelb_Thylacine_Refassem.4out -num_threads 64 -outfmt 4 -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $hybrid_assembly -out outfmt4/blastn."$query".GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly.4out -num_threads 64 -outfmt 4 -query $query

#OUT fmt 7
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $SRA_db -out outfmt7/blastn."$query".SRR5055303-6.blastDB.SRA.mview.out -num_threads 64 -outfmt "7 std qseq sseq stitle" -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $Refassem -out outfmt7/blastn."$query".GCA_007646695.1_UniMelb_Thylacine_Refassem.mview.out -num_threads 64 -outfmt "7 std qseq sseq stitle" -query $query
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db $hybrid_assembly -out outfmt7/blastn."$query".GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly.mview.out -num_threads 64 -outfmt "7 std qseq sseq stitle" -query $query

#MVIEW
mview -in blast -moltype dna -out fasta outfmt7/blastn."$query".SRR5055303-6.blastDB.SRA.mview.out|sed 's/^.*SRR/>SRR/g;s/[0-9]://g'|awk -F' ' '{print $1}'> mview/blastn."$query".SRR5055303-6.blastDB.SRA.mview.out.fasta
mview -in fasta -html head -find 'TGA:TAA:TAG' -css on  mview/blastn."$query".SRR5055303-6.blastDB.SRA.mview.out.fasta > mview/blastn."$query".SRR5055303-6.blastDB.SRA.mview.out.fasta.html

