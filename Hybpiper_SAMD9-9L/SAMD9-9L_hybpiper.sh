cat *_1.fq > XM_003771945.4_XM_031940019.1-SRR5055303-6_1.fq
cat *_2.fq > XM_003771945.4_XM_031940019.1-SRR5055303-6_2.fq

##blast
hybpiper assemble -r XM_003771945.4_XM_031940019.1-SRR5055303-6_1.fq XM_003771945.4_XM_031940019.1-SRR5055303-6_2.fq -t_dna SAMD9-SAMD9L.fa --prefix XM_003771945.4_XM_031940019.1-SRR5055303-6_blast --evalue 0.01 --cpu 64
#Without_thylacine
hybpiper assemble -r XM_003771945.4_XM_031940019.1-SRR5055303-6_1.fq XM_003771945.4_XM_031940019.1-SRR5055303-6_2.fq -t_dna SAMD9-SAMD9L_wo_thylacine.fa --prefix Without_thylacine_XM_003771945.4_XM_031940019.1-SRR5055303-6_blast --evalue 0.01 --cpu 64

hybpiper retrieve_sequences  --targetfile_dna SAMD9-SAMD9L.fa --single_sample_name XM_003771945.4_XM_031940019.1-SRR5055303-6_blast --fasta_dir XM_003771945.4_XM_031940019.1-SRR5055303-6_blast.extracted_seq dna
hybpiper retrieve_sequences  --targetfile_dna SAMD9-SAMD9L_wo_thylacine.fa --single_sample_name Without_thylacine_XM_003771945.4_XM_031940019.1-SRR5055303-6_blast --fasta_dir Without_thylacine_XM_003771945.4_XM_031940019.1-SRR5055303-6_blast.extracted_seq dna
