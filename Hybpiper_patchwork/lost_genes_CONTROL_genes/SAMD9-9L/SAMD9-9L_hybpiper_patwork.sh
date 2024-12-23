for i in *.fa
do
gene=`echo $i|sed 's/\.fa//g'`
#identify reads
blastn -task blastn  -evalue 0.05 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa -out blastn."$gene".1out -num_threads 64 -outfmt 1 -query $i -max_target_seqs 10000
#get read ids list
grep "^SRR" blastn."$gene".1out|grep "^SRR5055303"|grep "/1"|cut -f1,2 -d' ' |sort -u > "$gene".SRR5055303_1.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055303"|grep "/2"|cut -f1,2 -d' ' |sort -u > "$gene".SRR5055303_2.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055304"|grep "/1"|cut -f1,2 -d' ' |sort -u > "$gene".SRR5055304_1.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055304"|grep "/2"|cut -f1,2 -d' ' |sort -u > "$gene".SRR5055304_2.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055305"|cut -f1 -d' ' |sort -u > "$gene".SRR5055305.lst
grep "^SRR" blastn."$gene".1out|grep "^SRR5055306"|cut -f1 -d' ' |sort -u > "$gene".SRR5055306.lst
done

comm -23 XM_003771945.4.SRR5055303_1.lst XM_031940019.1.SRR5055303_1.lst|sort -u > XM_003771945.4.SRR5055303_1.uniq.lst
comm -23 XM_003771945.4.SRR5055303_2.lst XM_031940019.1.SRR5055303_2.lst|sort -u > XM_003771945.4.SRR5055303_2.uniq.lst
comm -23 XM_003771945.4.SRR5055304_1.lst XM_031940019.1.SRR5055304_1.lst|sort -u > XM_003771945.4.SRR5055304_1.uniq.lst
comm -23 XM_003771945.4.SRR5055304_2.lst XM_031940019.1.SRR5055304_2.lst|sort -u > XM_003771945.4.SRR5055304_2.uniq.lst
comm -23 XM_003771945.4.SRR5055305.lst XM_031940019.1.SRR5055305.lst|sort -u > XM_003771945.4.SRR5055305.uniq.lst
comm -23 XM_003771945.4.SRR5055306.lst XM_031940019.1.SRR5055306.lst|sort -u > XM_003771945.4.SRR5055306.uniq.lst

#extract reads
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303_1.fastq.gz XM_003771945.4.SRR5055303_1.uniq.lst > XM_003771945.4.SRR5055303_1.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303_2.fastq.gz XM_003771945.4.SRR5055303_2.uniq.lst > XM_003771945.4.SRR5055303_2.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055304_1.fastq.gz XM_003771945.4.SRR5055304_1.uniq.lst > XM_003771945.4.SRR5055304_1.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055304_2.fastq.gz XM_003771945.4.SRR5055304_2.uniq.lst > XM_003771945.4.SRR5055304_2.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055305_1.fastq.gz XM_003771945.4.SRR5055305.uniq.lst > XM_003771945.4.SRR5055305_1.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055305_2.fastq.gz XM_003771945.4.SRR5055305.uniq.lst > XM_003771945.4.SRR5055305_2.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055306_1.fastq.gz XM_003771945.4.SRR5055306.uniq.lst > XM_003771945.4.SRR5055306_1.uniq.fq
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055306_2.fastq.gz XM_003771945.4.SRR5055306.uniq.lst > XM_003771945.4.SRR5055306_2.uniq.fq

#hybpiper

hybpiper assemble -r XM_003771945.4.SRR5055303_1.uniq.fq XM_003771945.4.SRR5055303_2.uniq.fq -t_dna XM_003771945.4.fa --prefix XM_003771945.4.SRR5055303_blast --evalue 0.05 --cpu 64
hybpiper assemble -r XM_003771945.4.SRR5055304_1.uniq.fq XM_003771945.4.SRR5055304_2.uniq.fq -t_dna XM_003771945.4.fa --prefix XM_003771945.4.SRR5055304_blast --evalue 0.05 --cpu 64
hybpiper assemble -r XM_003771945.4.SRR5055305_1.uniq.fq XM_003771945.4.SRR5055305_2.uniq.fq -t_dna XM_003771945.4.fa --prefix XM_003771945.4.SRR5055305_blast --evalue 0.05 --cpu 64
hybpiper assemble -r XM_003771945.4.SRR5055306_1.uniq.fq XM_003771945.4.SRR5055306_2.uniq.fq -t_dna XM_003771945.4.fa --prefix XM_003771945.4.SRR5055306_blast --evalue 0.05 --cpu 64

# retrive sequence

hybpiper retrieve_sequences  --targetfile_dna XM_003771945.4.fa --single_sample_name XM_003771945.4.SRR5055303_blast --fasta_dir SRR5055303.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna XM_003771945.4.fa --single_sample_name XM_003771945.4.SRR5055304_blast --fasta_dir SRR5055304.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna XM_003771945.4.fa --single_sample_name XM_003771945.4.SRR5055305_blast --fasta_dir SRR5055305.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna XM_003771945.4.fa --single_sample_name XM_003771945.4.SRR5055306_blast --fasta_dir SRR5055306.extracted_seq dna

#patchwork
cat XM_003771945.4.SRR5055303_1.uniq.lst XM_003771945.4.SRR5055303_2.uniq.lst XM_003771945.4.SRR5055304_1.uniq.lst XM_003771945.4.SRR5055304_2.uniq.lst XM_003771945.4.SRR5055305.uniq.lst XM_003771945.4.SRR5055306.uniq.lst > XM_003771945.4.uniq.lst
faTrans -stop XM_003771945.4.fa XM_003771945.4.faTrans.fa
sed -i 's/>Sarcophilus_harrisii-/>Sarcophilus_harrisii@/g' XM_003771945.4.faTrans.fa
seqtk subseq /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/SRA/SRR5055303-6.blastDB.fa XM_003771945.4.uniq.lst > XM_003771945.patchwork.fasta
julia --project=/media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/tools/patchwork /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Hybpiper/tools/patchwork/src/Patchwork.jl --contigs XM_003771945.patchwork.fasta --reference XM_003771945.4.faTrans.fa --output-dir XM_003771945.4_patchwork

