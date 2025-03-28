clade=$1
cd $clade
mkdir RELAX
cp "$clade".aln "$clade".nwk RELAX
cd RELAX
tree="$clade".nwk
for i in "$clade".aln
do
t=`grep ">" $i|wc -l`
grep ">" $i|sed 's/>//g' > taxlist.txt
for j in `cat taxlist.txt`
do
sed "s/$j/$j{fg}/g" $tree > "$j"_treeLabled.txt
echo -e  "1\n7\n1\n"$PWD"/$i\n"$PWD"/"$j"_treeLabled.txt\n2\n2" > "$j"_tree.config
/media/morpheus/sagar/BUDDHA/HYPHY_PAML_gBGC/hyphy-2.5.43rc/HYPHYMP <"$j"_tree.config > "$j"_treeoutput_relax
done
done


### taking HYPHYMP output result in one file
echo -e "test back pval kval"|sed 's/ /\t/g' > HYPHY_RELAX.Results.txt
for d in `ls -1 *_relax`
do
pval=`grep "^Like" $d|awk '{print $6}'|sed 's/\*\*\.//g'|head -1`
kval1=`grep "Relaxation/intensification" $d|awk '{print $6}'|head -1`
test=`grep "_Test_ set:" $d|awk '{print $9}'|sed 's/\`//g'`
back=`grep "_Reference_ set:" $d|cut -f2 -d ":"|sed 's/ //g'`

echo -e "$test $back $pval $kval1"|sed 's/ /\t/g' >> HYPHY_RELAX.Results.txt
done
cut -f1,3,4 HYPHY_RELAX.Results.txt > HYPHY_RELAX.Results.sorted.txt
cd ../..


