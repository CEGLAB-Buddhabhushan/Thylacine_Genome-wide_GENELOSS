###Monodelphis_domestica
klumpy scan_alignments --alignment_map Monodelphis_domestica.merge.sorted.bam --threads 16
klumpy alignment_plot --alignment_map Monodelphis_domestica.merge.sorted.bam --reference NC_077232.1 --candidates Monodelphis_domestica.merge.sorted_Candidate_Regions.tsv --min_len 3000 --window_size 10000 --window_step 5000 --color red --vertical_line_gaps --vertical_line_klumps --format svg --leftbound 266519535 --rightbound 266794692 --annotation GCF_027887165.1_mMonDom1.pri_genomic.gtf --gap_file Monodelphis_domestica.HSD17B13.NC_077232.1_gaps.tsv
#Found the following genes: IBSP NUDT9 LOC103092505 LOC100016979 KLHL8
