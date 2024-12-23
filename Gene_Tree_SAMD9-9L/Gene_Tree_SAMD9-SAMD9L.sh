mkdir Gene_Tree_SAMD9-9L
cd Gene_Tree_SAMD9-9L
sed 's/>/>SAMD9_/g' ../SAMD9/SAMD9.fa > SAMD9-SAMD9L.fa
sed 's/>/>SAMD9L_/g' ../SAMD9L/SAMD9L.fa >> SAMD9-SAMD9L.fa
sed -i 's/-//g' SAMD9-SAMD9L.fa
guidance=/home/morpheus/gprc6a/guidance.v2.02/www/Guidance/guidance.pl
perl $guidance --program GUIDANCE --seqFile SAMD9-SAMD9L.fa --msaProgram PRANK --seqType nuc --outDir SAMD9-SAMD9L.100_PRANK --genCode 1 --bootstraps 100 --proc_num 48
seqkit sort SAMD9-SAMD9L.100_PRANK/MSA.PRANK.aln.With_Names >SAMD9-SAMD9L.PRANK.aln
rm -r SAMD9-SAMD9L.100_PRANK

mkdir iqtree2_MF
cp SAMD9-SAMD9L.PRANK.aln iqtree2_MF
cd iqtree2_MF
#One would be more confident if a clade has its SH-aLRT >= 80% and UFboot >= 95%.
/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.PRANK.aln -m MFP --alrt 1000 -B 1000 --boot-trees


##phylogenetic tree for 1kb region and Gene concordance factor (gCF)
mkdir splitMSA_1kb
cp SAMD9-SAMD9L.PRANK.aln splitMSA_1kb
cd splitMSA_1kb
seqkit sliding -s 1000 -W 1000 -o output_chunks.fasta SAMD9-SAMD9L.PRANK.aln
/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.PRANK.aln -m MFP --alrt 1000 -B 1000 --boot-trees
grep "SAMD9L_Dasyurus_viverrinus_sliding" output_chunks.fasta|cut -f2 -d':' > window.lst
seqtk seq output_chunks.fasta > SAMD9-SAMD9L.1kb_window.aln
for i in `cat window.lst`
do
grep "$i" SAMD9-SAMD9L.1kb_window.aln |sed 's/>//g' > "$i".lst
seqtk subseq SAMD9-SAMD9L.1kb_window.aln "$i".lst  > SAMD9-SAMD9L.1kb_window."$i".aln
sed -i 's/_sliding:.*//' SAMD9-SAMD9L.1kb_window."$i".aln
/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.1kb_window."$i".aln -m MFP --alrt 1000 -B 1000 --boot-trees
#Gene concordance factor (gCF)
/home/morpheus/anaconda3/bin/iqtree2 -t SAMD9-SAMD9L.PRANK.aln.treefile --gcf SAMD9-SAMD9L.1kb_window."$i".aln.treefile --prefix gCF."$i".concord
done
