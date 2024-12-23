export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export NXF_VER=22.10.0

chr=$1
cd $chr
for i in "$chr"_chr.fa
do
target_name=`echo Tasmanian_devil_"$chr"`
faToTwoBit "$chr"_chr.fa "$chr"_chr.2bit
/media/morpheus/sagar/BUDDHA/TOGA/toga.py ./Chain_"$target_name"_Tasmanian_wolf/"$target_name".Tasmanian_wolf.final.chain.gz GCF_902635505.1_mSarHar1.11_genomic."$chr"_chr.bed12 "$chr"_chr.2bit ../Tasmanian_wolf.2bit --kt --pn TOGA_"$target_name"_Tasmanian_wolf_o2o --nc /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/ --cb 3,5,15,25,50 --cjn 500 --nb /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/cesar_bigmem_config.nf --ces --o2o --chain_jobs_num 60
done
cd ..

