job_no=$1
cd "$job_no"
for i in *.fa
do
gene=`echo $i|sed 's/\.fa//g'`
#identify reads
blastn -task blastn  -evalue 0.05 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out blastn."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^SRR" blastn."$gene".1out|grep "^SRR5055303"|grep "/1"|cut -f1,2 -d' '  > SRR5055303_1.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055303"|grep "/2"|cut -f1,2 -d' '  > SRR5055303_2.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055304"|grep "/1"|cut -f1,2 -d' '  > SRR5055304_1.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055304"|grep "/2"|cut -f1,2 -d' '  > SRR5055304_2.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055305"|cut -f1 -d' '  > SRR5055305.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055306"|cut -f1 -d' '  > SRR5055306.lst

#extract reads
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303_1.fastq.gz SRR5055303_1.lst > "$gene".SRR5055303_1.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303_2.fastq.gz SRR5055303_2.lst > "$gene".SRR5055303_2.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055304_1.fastq.gz SRR5055304_1.lst > "$gene".SRR5055304_1.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055304_2.fastq.gz SRR5055304_2.lst > "$gene".SRR5055304_2.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055305_1.fastq.gz SRR5055305.lst > "$gene".SRR5055305_1.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055305_2.fastq.gz SRR5055305.lst > "$gene".SRR5055305_2.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055306_1.fastq.gz SRR5055306.lst > "$gene".SRR5055306_1.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055306_2.fastq.gz SRR5055306.lst > "$gene".SRR5055306_2.fq

#hybpiper

hybpiper assemble -r "$gene".SRR5055303_1.fq "$gene".SRR5055303_2.fq -t_dna $i --prefix "$gene".SRR5055303_blast --evalue 0.05 --cpu 64
hybpiper assemble -r "$gene".SRR5055304_1.fq "$gene".SRR5055304_2.fq -t_dna $i --prefix "$gene".SRR5055304_blast --evalue 0.05 --cpu 64
hybpiper assemble -r "$gene".SRR5055305_1.fq "$gene".SRR5055305_2.fq -t_dna $i --prefix "$gene".SRR5055305_blast --evalue 0.05 --cpu 64
hybpiper assemble -r "$gene".SRR5055306_1.fq "$gene".SRR5055306_2.fq -t_dna $i --prefix "$gene".SRR5055306_blast --evalue 0.05 --cpu 64

# retrive sequence

hybpiper retrieve_sequences  --targetfile_dna $i --single_sample_name "$gene".SRR5055303_blast --fasta_dir SRR5055303.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna $i --single_sample_name "$gene".SRR5055304_blast --fasta_dir SRR5055304.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna $i --single_sample_name "$gene".SRR5055305_blast --fasta_dir SRR5055305.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna $i --single_sample_name "$gene".SRR5055306_blast --fasta_dir SRR5055306.extracted_seq dna

#patchwork

faTrans -stop $i "$gene".faTrans.fa
sed -i 's/>Sarcophilus_harrisii-/>Sarcophilus_harrisii@/g' "$gene".faTrans.fa
grep "^SRR" blastn."$gene".1out | cut -f1 -d' ' > "$gene".patchwork.lst
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa "$gene".patchwork.lst > "$gene".patchwork.fasta
julia --project=/media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/tools/patchwork /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/tools/patchwork/src/Patchwork.jl --contigs "$gene".patchwork.fasta --reference "$gene".faTrans.fa --output-dir "$gene"_patchwork

echo "HYBPIPER and PATCHWORK completed for $gene"
done
cd ..
