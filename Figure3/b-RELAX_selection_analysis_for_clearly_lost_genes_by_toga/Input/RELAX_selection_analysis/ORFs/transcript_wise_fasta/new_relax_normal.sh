#!/bin/bash
#conda activate hyphy_env
clade=$1
cd $clade
tree="$clade".nwk
sed 's/Thylacinus_cynocephalus/Thylacinus_cynocephalus{fg}/g' $tree > Thylacinus_cynocephalus_treeLabled.txt
hyphy relax --alignment "$clade".aln --tree Thylacinus_cynocephalus_treeLabled.txt --test fg --output Thylacinus_cynocephalus_treeoutput_relax
