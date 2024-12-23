#create a folder with gene name and save the cds and exons in that folder with suffix .fa
for i in *.fa; do  ./../BLASTn_fmt1_4_6_7_mview.sh $i; done
