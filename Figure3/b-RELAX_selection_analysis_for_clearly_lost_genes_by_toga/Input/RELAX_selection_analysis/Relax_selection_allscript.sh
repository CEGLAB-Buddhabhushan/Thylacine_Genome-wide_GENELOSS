for i in `cat Gene_status_not_Devil-specific.tsv`
do
cd $i
for toga in TOGA_*
do
ref=`echo $toga|cut -f2,3 -d'_'`
query=`echo $toga|cut -f4,5 -d'_'`
cat $toga/codon.fasta |sed 's/ //g;s/X/N/g;s/-//g'|sed "s/>.*REFERENCE/>$ref/g;s/>.*QUERY/>$query/g" >> /media/morpheus/disk1/BUDDHA_merged_data/history_of_gene_loss/Sh_L_transcript/RELAX_selection_analysis/ORFs/"$i".codon.fasta
done
cd ..
done
cd RELAX_selection_analysis/ORFs/
for i in `cut -f1 ../../Gene_status_not_Devil-specific.tsv`
do
sed -e ':a' -e 'N' -e '/\n\s*$/!P; D; ba'  "$i".codon.fasta >  "$i".codon.reformatted.fasta
mafft "$i".codon.reformatted.fasta > "$i".codon.reformatted.mafft.fasta
mview -in fasta -out fasta -reference 1 -moltype dna -pcid reference -sort cov  -minident 70 "$i".codon.reformatted.mafft.fasta|cdskit aggregate |cut -f1 -d' '|seqtk seq|sed 's/^---/atg/' |sed 's/-//g'|grep -B1 '^atg.*nnn$'|sed 's/nnn$//'|awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' > "$i".codon.reformatted.mafft.mview.fasta
perl multifasta.pl "$i".codon.reformatted.mafft.mview.fasta > "$i".codon.reformatted.mafft.mview.multifasta.fasta
guidance=/home/morpheus/gprc6a/guidance.v2.02/www/Guidance/guidance.pl
perl $guidance --program GUIDANCE --seqFile "$i".codon.reformatted.mafft.mview.multifasta.fasta --msaProgram PRANK --seqType codon --outDir "$i".100_PRANK --genCode 1 --bootstraps 100 --proc_num 48
seqkit sort "$i".100_PRANK/MSA.PRANK.aln.With_Names >"$i"_PRANK.aln
rm -r "$i".100_PRANK
done


cat *.aln|grep ">"|sed 's/>//g'|sort -u > Species.lst
#get the species tree from time tree website
#((((Macrotis_lagotis:58.00000000,(Thylacinus_cynocephalus:38.47000000,(Myrmecobius_fasciatus:33.18550000,(Sminthopsis_crassicaudata:24.59895000,(Antechinus_flavipes:18.03550000,(Sarcophilus_harrisii:9.44000000,Dasyurus_viverrinus:9.44000000)'14':8.59550000)'13':6.56345000)'25':8.58655000)'37':5.28450000)'36':19.53000000)'35':2.77164000,((((Phalanger_gymnotis:25.35000000,Trichosurus_vulpecula:25.35000000)'34':19.95000000,(Potorous_gilbertii:21.31000000,((Notamacropus_eugenii:7.37000000,(Macropus_fuliginosus:2.11000000,Macropus_giganteus:2.11000000)'43':5.26000000)'33':0.95500000,Lagorchestes_hirsutus:8.32500000)'51':12.98500000)'50':23.99000000)'49':2.75500000,((Pseudochirops_cupreus:6.50500000,Pseudochirops_corinnae:6.50500000)'57':16.90000000,Pseudocheirus_occidentalis:23.40500000)'60':24.65000000)'56':4.58500000,Vombatus_ursinus:52.64000000)'48':8.13164000)'63':2.22836000,Dromiciops_gliroides:63.00000000)'47':15.18680000,(Monodelphis_domestica:29.14096000,Gracilinanus_agilis:29.14096000)'68':49.04584000);


#### remove internodes and unroot tree
for i in Species.nwk
do
sed -e "s/'[^()]*'//g" $i > Species.no_internodes.nwk
echo 'library(ape)' > tree_script.r
echo "a<-read.tree(\"Species.no_internodes.nwk\")" >> tree_script.r
echo 'b<-unroot(a)' >> tree_script.r
echo "write.tree(b,\"Species.no_internodes.unroot.nwk\")" >> tree_script.r
Rscript tree_script.r
done
# delete the empty .aln files
find . -name "*.aln" -size 0 -exec rm {} \;
##keep tips
for aln in *_PRANK.aln
do
transcriptID=`echo $aln |sed 's/_PRANK.aln//g'`
tips=`grep ">" "$aln" | sed 's/>//g' | awk 'BEGIN { printf("tip <- c(\n") } { printf("    \"%s\", ", $0) } END { printf("\b\b\n)\n") }'|sed 's/ //g'|tr '\n' ' '|sed 's/ //g;s/)/")/g'` 
echo 'library(ape)' > keep_tip.r
echo 'a <- read.tree("Species.no_internodes.nwk")' >> keep_tip.r
echo "$tips" >> keep_tip.r
echo 'b <- keep.tip(a, tip)' >> keep_tip.r
echo 'c<-unroot(b)' >> keep_tip.r
echo 'write.tree(c, "transcriptID.tree")' >> keep_tip.r
sed -i 's/\x08//g;s/,")/)/g' keep_tip.r
Rscript keep_tip.r
mv transcriptID.tree "$transcriptID".nwk 
echo $transcriptID
done

mkdir transcript_wise_fasta
for aln in *_PRANK.aln
do
transcriptID=`echo $aln |sed 's/_PRANK.aln//g'`
mkdir -p transcript_wise_fasta/$transcriptID
mv $aln transcript_wise_fasta/$transcriptID/"$transcriptID".aln
mv "$transcriptID".nwk transcript_wise_fasta/$transcriptID/"$transcriptID".nwk
done
cd transcript_wise_fasta/

for dir in `ls -d */`
do
transcriptID=`echo $dir |sed 's/[/]//g'`
thylacine="no"
    if grep -q "Thylacinus_cynocephalus" $transcriptID/"$transcriptID".aln; then
        thylacine="yes"
    fi
echo -e "$transcriptID\t$thylacine" >> Thylacine_status.list
done

#add lost gene directory with aln and tree; add this genes in Thylacine_status.list
#run relax
#conda activate hyphy_env
for i in `grep "yes" Thylacine_status.list|cut  -f1 `
do
./new_relax_normal_wrapper.sh $i
done


### taking HYPHYMP output result in one file
echo -e "Gene\tTest\tLRT\tp-value\tK" > HYPHY_RELAX.Results.txt
for i in `grep "yes" Thylacine_status.list|cut  -f1 `
do
cd $i
d="Thylacinus_cynocephalus_treeoutput_relax"
test=`grep "Test" $d |tail -1|cut -f1 -d':'|sed 's/"//g;s/ //g'`
lrt=`grep "LRT" $d |cut -f2 -d':'|sed 's/,//g'`
pval=`grep "p-value" $d |cut -f2 -d':'|sed 's/,//g'`
kval=`grep "relaxation or intensification parameter" $d |cut -f2 -d':'|sed 's/,//g'`
echo -e "$i\t$test\t$lrt\t$pval\t$kval" >> ../HYPHY_RELAX.Results.txt
cd ..
done
awk 'NF != 5 {print $1}' HYPHY_RELAX.Results.txt  > hyphy_error.lst
for i in `cat hyphy_error.lst`
do
bash new_relax_TOLERATE_NUMERICAL_ERRORS.sh $i
done
echo -e "Gene\tTest\tLRT\tp-value\tK" > HYPHY_RELAX.Results.txt
for i in `grep "yes" Thylacine_status.list|cut  -f1 `
do
cd $i
d="Thylacinus_cynocephalus_treeoutput_relax"
test=`grep "Test" $d |tail -1|cut -f1 -d':'|sed 's/"//g;s/ //g'`
lrt=`grep "LRT" $d |cut -f2 -d':'|sed 's/,//g'`
pval=`grep "p-value" $d |cut -f2 -d':'|sed 's/,//g'`
kval=`grep "relaxation or intensification parameter" $d |cut -f2 -d':'|sed 's/,//g'`
echo -e "$i\t$test\t$lrt\t$pval\t$kval" >> ../HYPHY_RELAX.Results.txt
cd ..
done













