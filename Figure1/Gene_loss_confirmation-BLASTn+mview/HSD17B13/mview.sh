#gt -> ga
mview -in fasta -html head -css on -range 195:228 -hide all  -show '1,3,2,18,24,34,45,46'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_HSD17B13_XM_003774508.4_Exon1_intron1_25nt.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.fasta > HSD17B13_gt-ga.sorted.html

#-1 (G)
mview -in fasta -html head -css on -range 562:590 -hide all  -show '1,2,3,23,29,35,36,499,318' blastn.Sarcophilus_harrisii_HSD17B13_XM_003774508.4.CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.fa.mview.fasta > HSD17B13_-1_codon193.sorted.html
