##generate input for selection analysis
gene=$1
# to remove internode labels from the TimeTree nwk file
for i in "$gene".nwk
do
sed -e "s/'[^()]*'//g" $i > temp.nwk
mv temp.nwk $i
echo 'library(ape)' > tree_script.r
echo "a<-read.tree(\"$i\")" >> tree_script.r
echo 'b<-unroot(a)' >> tree_script.r
echo "write.tree(b,\"$i.tree\")" >> tree_script.r
Rscript tree_script.r
mv $i.tree $i
done

# This script check whether the species name in the fasta file and nwk file are same or not.
for i in "$gene".ORF.fa
do
grep ">" $i|sed 's/>//g' > $i.txt
j="$gene".nwk
sed 's/(/\n/g' $j|sed 's/)/\n/g' |sed 's/;/\n/g' |sed 's/:/\n/g' |sed 's/,/\n/g' |grep "^[A-Z]" > $j.txt
echo $j
cat $i.txt $j.txt |sort|uniq -c |awk '$1<2 {print $2}'
rm *.txt
done

guidance=/home/ceglab358/BUDDHA/Tools/guidance.v2.02/www/Guidance/guidance.pl
for i in `ls "$gene".ORF.fa`
do
j=`echo $i|sed 's/.ORF.fa//g'`
perl $guidance --program GUIDANCE --seqFile "$i" --msaProgram PRANK --seqType codon --outDir "$i".100_PRANK --genCode 1 --bootstraps 100 --proc_num 16
cp "$i".100_PRANK/MSA.PRANK.aln.With_Names "$j".aln
rm -r "$i".100_PRANK
done

############## check gene tree of SAMD9 and SAMD9L
mkdir SAMD9-SAMD9L
cd SAMD9-SAMD9L
sed 's/>/>SAMD9_/g' ../SAMD9/SAMD9.fa > SAMD9-SAMD9L.fa
sed 's/>/>SAMD9L_/g' ../SAMD9L/SAMD9L.fa >> SAMD9-SAMD9L.fa
sed -i 's/-//g' SAMD9-SAMD9L.fa
guidance=/home/morpheus/gprc6a/guidance.v2.02/www/Guidance/guidance.pl
perl $guidance --program GUIDANCE --seqFile SAMD9-SAMD9L.fa --msaProgram PRANK --seqType nuc --outDir SAMD9-SAMD9L.100_PRANK --genCode 1 --bootstraps 100 --proc_num 48
perl $guidance --program GUIDANCE --seqFile SAMD9-SAMD9L.fa --msaProgram MUSCLE --seqType nuc --outDir SAMD9-SAMD9L.100_MUSCLE --genCode 1 --bootstraps 100 --proc_num 48
perl $guidance --program GUIDANCE --seqFile SAMD9-SAMD9L.fa --msaProgram CLUSTALW --seqType nuc --outDir SAMD9-SAMD9L.100_CLUSTALW --genCode 1 --bootstraps 100 --proc_num 48
perl $guidance --program GUIDANCE --seqFile SAMD9-SAMD9L.fa --msaProgram MAFFT --seqType nuc --outDir SAMD9-SAMD9L.100_MAFFT --genCode 1 --bootstraps 100 --proc_num 48
seqkit sort SAMD9-SAMD9L.100_PRANK/MSA.PRANK.aln.With_Names >SAMD9-SAMD9L.PRANK.aln
seqkit sort SAMD9-SAMD9L.100_MUSCLE/MSA.MUSCLE.aln.With_Names >SAMD9-SAMD9L.MUSCLE.aln
seqkit sort SAMD9-SAMD9L.100_CLUSTALW/MSA.CLUSTALW.aln.With_Names >SAMD9-SAMD9L.CLUSTALW.aln
seqkit sort SAMD9-SAMD9L.100_MAFFT/MSA.MAFFT.aln.With_Names >SAMD9-SAMD9L.MAFFT.aln
/media/morpheus/sagar/BUDDHA/Tools/mumsa-1.0/mumsa -g -s SAMD9-SAMD9L.PRANK.aln SAMD9-SAMD9L.MUSCLE.aln SAMD9-SAMD9L.CLUSTALW.aln SAMD9-SAMD9L.MAFFT.aln >> SAMD9-SAMD9L.log
/media/morpheus/sagar/BUDDHA/Tools/mumsa-1.0/mumsa SAMD9-SAMD9L.PRANK.aln SAMD9-SAMD9L.MUSCLE.aln SAMD9-SAMD9L.CLUSTALW.aln SAMD9-SAMD9L.MAFFT.aln -a -q > SAMD9-SAMD9L.highest_MOS.NT.aln
rm -r SAMD9-SAMD9L.100_PRANK SAMD9-SAMD9L.100_MUSCLE SAMD9-SAMD9L.100_CLUSTALW SAMD9-SAMD9L.100_MAFFT
mkdir raxml-ng 
cp SAMD9-SAMD9L.highest_MOS.NT.aln raxml-ng
mkdir iqtree2
cp SAMD9-SAMD9L.highest_MOS.NT.aln iqtree2
modeltest-ng -i SAMD9-SAMD9L.highest_MOS.NT.aln -o SAMD9-SAMD9L.highest_MOS.NT.aln.modeltest
model=$(awk '{if($1=="BIC"){print $2}}' *.modeltest.log| tail -1)
echo $model
raxml-ng --all --msa SAMD9-SAMD9L.highest_MOS.NT.aln --model $model --bs-trees 1000 --threads 48 --prefix boot 
#sed -e "s/'[^()]*'//g" boot.raxml.bestTree > SAMD9-SAMD9L.raxml.bestTree.nwk
raxmlHPC -f b -m PROTGAMMAILG -n output_bootstrap.tree -t boot.raxml.bestTree -z boot.raxml.bootstraps


/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.highest_MOS.NT.aln -st DNA --boot 1000
sed -e "s/'[^()]*'//g" "$gene".highest_MOS.NT.aln.contree > "$gene".iqtree.contree.nwk
#Consensus construction and bootstrap value assignment
iqtree -con mytrees -minsup 0.5
iqtree -net mytrees

mkdir iqtree2_MF
cp SAMD9-SAMD9L.highest_MOS.NT.aln iqtree2_MF
cd iqtree2_MF
#One would be more confident if a clade has its SH-aLRT >= 80% and UFboot >= 95%.
/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.highest_MOS.NT.aln -m MFP --alrt 1000 -B 1000 --boot-trees
#Gene concordance factor (gCF)
/home/morpheus/anaconda3/bin/iqtree2 -t concat.treefile --gcf loci.treefile --prefix concord
#Site concordance factor (sCF)
/home/morpheus/anaconda3/bin/iqtree2 -te concat.treefile -s ALN_FILE --scfl 100 --prefix concord


##phylogenetic tree for 1kb region
mkdir splitMSA_1kb
cp SAMD9-SAMD9L.highest_MOS.NT.aln splitMSA_1kb
cd splitMSA_1kb
seqkit sliding -s 1000 -W 1000 -o output_chunks.fasta SAMD9-SAMD9L.highest_MOS.NT.aln
/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.highest_MOS.NT.aln -m MFP --alrt 1000 -B 1000 --boot-trees
grep "SAMD9L_Dasyurus_viverrinus_sliding" output_chunks.fasta|cut -f2 -d':' > window.lst
seqtk seq output_chunks.fasta > SAMD9-SAMD9L.highest_MOS.NT.1kb.aln
for i in `cat window.lst`
do
grep "$i" SAMD9-SAMD9L.highest_MOS.NT.1kb.aln |sed 's/>//g' > "$i".lst
seqtk subseq SAMD9-SAMD9L.highest_MOS.NT.1kb.aln "$i".lst  > SAMD9-SAMD9L.highest_MOS.NT."$i".aln
sed -i 's/_sliding:.*//' SAMD9-SAMD9L.highest_MOS.NT."$i".aln
/home/morpheus/anaconda3/bin/iqtree2 -T AUTO -s SAMD9-SAMD9L.highest_MOS.NT."$i".aln -m MFP --alrt 1000 -B 1000 --boot-trees
#Gene concordance factor (gCF)
/home/morpheus/anaconda3/bin/iqtree2 -t SAMD9-SAMD9L.highest_MOS.NT.aln.treefile --gcf SAMD9-SAMD9L.highest_MOS.NT."$i".aln.treefile --prefix gCF."$i".concord
done
mkdir align
cp SAMD9-SAMD9L.highest_MOS.NT.1-1000.aln SAMD9-SAMD9L.highest_MOS.NT.1001-2000.aln SAMD9-SAMD9L.highest_MOS.NT.2001-3000.aln SAMD9-SAMD9L.highest_MOS.NT.3001-4000.aln align/
#Site concordance factor (sCF)
/home/morpheus/anaconda3/bin/iqtree2 -t SAMD9-SAMD9L.highest_MOS.NT.aln.treefile -p align/ --scfl 100 --prefix sCFconcord

