#MVIEW
for i in SAMD9L
do
cd $i
query="Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa"
mview -in blast -moltype dna -out fasta outfmt7/blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_UniMelb_Thylacine_Refassem.mview.out|sed 's/^.*VAHE/>VAHE/g;s/[0-9]://g'|awk -F' ' '{print $1}' > main_fig/blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta
mview -in blast -moltype dna -out fasta outfmt7/blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly.mview.out|sed 's/^.*CM/>CM/g;s/[0-9]://g'|awk -F' ' '{print $1}' >> main_fig/blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta
mview -in blast -moltype dna -out fasta outfmt7/blastn."$query".SRR5055303-6.blastDB.SRA.mview.out|sed 's/^.*SRR/>SRR/g;s/[0-9]://g'|awk -F' ' '{print $1}' >> main_fig/blastn."$query".GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta
cd ..
done
cd SAMD9L/main_fig
mview -in fasta -html head -find 'TGA:TAA:TAG' -css on  blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta.html


#subset the location and reads
#-1
mview -in fasta -html head -css on -range 840:875 -hide all  -show '1,2,4,280,548,636,646,707,1151,1167'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_-1.SRA.mview.sorted.html
#CGA> TGA
mview -in fasta -html head -css on -range 1520:1540 -hide all  -show '1,2,4,102,140,170,171,189,235,1551,1570,1572'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_CGA-TGA.SRA.mview.sorted.html
#-4
mview -in fasta -html head -css on -range 1925:1945 -hide all  -show '1,2,4,147,211,215,236,340,444,1402,1448,1420'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_-4_SRA.mview.sorted.html
#+1
mview -in fasta -html head -css on -range 1990:1205 -hide all  -show '1,2,4,58,77,78,95,114,148'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_+1_SRA.mview.sorted.html
#CGA> TAG
mview -in fasta -html head -css on -range 2025:2055 -hide all  -show '1,2,4,921,934,935,962,988'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_CGA-TAG_SRA.mview.sorted.html
#CAG> TAG
mview -in fasta -html head -css on -range 3610:3630 -hide all  -show '1,2,4,184,192,201,223,268,357,437'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_CAG-TAG_SRA.mview.sorted.html
#CAA> TAA
mview -in fasta -html head -css on -range 4375:4405 -hide all  -show '1,2,4,45,54,69,108,251,557,1597'  -label0 -label2 -label3 -label4 -label5 -label6 -label7 -label8 blastn.Sarcophilus_harrisii_LOC100924681_XM_031940019.1_CDS.fa.GCA_007646695.1_GCA_007646695.3_SRR5055303-6.blastDB.SRA.mview.out.fasta > SAMD9L_CAA-TAA_SRA.mview.sorted.html


