chr=$1
cd $chr
ulimit -n 16384
for i in "$chr"_chr.fa
do
sed -i 's/ .*//;s/\./_/' $i
target_name=`echo Tasmanian_devil_"$chr"`
query_name=Tasmanian_wolf
/media/morpheus/sagar/BUDDHA/TOGA_new/make_lastz_chains/make_chains.py $target_name $query_name $i ../GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly_genomic.fna --pd Chain_"$target_name"_"$query_name" -f --chaining_memory 30 
done
