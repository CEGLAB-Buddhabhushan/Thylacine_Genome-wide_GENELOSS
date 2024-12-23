for gene in SAMD9L VWA7 CUZD1 HSD17B13 19spCUZD1
do
cd $gene/Geneloss_timing

echo -e "Foreground\tBackgrounds\tmodel\tTm\twf\twm\twp\tTf\tTp\tTf2\tTp2" > "$gene".Gene_loss_timing.tsv
for ctl in `ls *ctl`
do
#codeml "$ctl"
tree=`ls -1 *nwk`
#sed -e 's/(/\n/g' -e 's/)/\n/g' -e 's/:/\n/g' -e 's/,/\n/g' $tree |grep "^[A-Z]" |grep "#" |sed -e 's/#.*//g'  -e 's/ //g' > fg.txt 
#sed -e 's/(/\n/g' -e 's/)/\n/g' -e 's/:/\n/g' -e 's/,/\n/g' $tree |grep "^[A-Z]" |grep -v "#"  > bg.txt 
outfile=`grep "outfile =" $ctl|awk '{print $3}'`
for sp in `cat fg.txt`
do
bg=`cat bg.txt|tr '\n' ','|sed 's/,$/\n/g'`
for div_time in 38.2 38.5
do
tm=`echo $div_time`
mdl=`grep "Codon frequency model:" $outfile|awk -F ":" '{print $2}'|awk '{print $1}'`
wf=`grep "w (dN/dS) for branches:" $outfile|awk '{print $5}'`
wm=`grep "w (dN/dS) for branches:" $outfile|awk '{print $6}'`
wp=`grep "w (dN/dS) for branches:" $outfile|awk '{print $7}'`
tf=`echo $tm $wm $wf | awk '{print $1*(($2 - 1)/($3 - 1))}'`
tp=`echo $tm $tf|awk '{print $1 - $2}'`
tf2=`echo $tm $tf $tp|awk '{print ($1*$2)/($2+(0.7*$3))}'`
tp2=`echo $tm $tf2|awk '{print $1-$2}'`
echo -e "$sp\t$bg\t$mdl\t$tm\t$wf\t$wm\t$wp\t$tf\t$tp\t$tf2\t$tp2" >> "$gene".Gene_loss_timing.tsv
done
#rm codeml.ctl
done
done
cp "$gene".Gene_loss_timing.tsv /media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/Thylacine_selection/Geneloss_timing
cd ../..
echo $gene
done



awk '{print $1,$3,$8,$10}' OFS='\t' outfilecompile|sed '1d'> outfilecompile_combined

Rscript Pseudogenization_time_compile.r

