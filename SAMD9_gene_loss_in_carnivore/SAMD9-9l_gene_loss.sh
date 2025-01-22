curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat

# Download human genome and gtf
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz && wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

## download datasets for HEPACAM2:253012 and CDK6:1021
for order in Carnivora Perissodactyla Artiodactyla
do
for gene in HEPACAM2_253012 CDK6_1021
do
gene_id=`echo $gene|cut -f2 -d'_'`
gene_name=`echo $gene|cut -f1 -d'_'`
./datasets download gene gene-id $gene_id --ortholog $order --no-progressbar --include product-report --filename "$gene_name"_"$order".zip
rm -rf "$gene_name"_"$order"
unzip -q -o "$gene_name"_"$order".zip -d "$gene_name"_"$order"
done
done

for order in Pholidota
do
for gene in HEPACAM2_253012 CDK6_1021
do
gene_id=`echo $gene|cut -f2 -d'_'`
gene_name=`echo $gene|cut -f1 -d'_'`
./datasets download gene gene-id $gene_id --ortholog 9971 --no-progressbar --include product-report --filename "$gene_name"_"$order".zip
rm -rf "$gene_name"_"$order"
unzip -q -o "$gene_name"_"$order".zip -d "$gene_name"_"$order"
done
done

for i in CDK6_* HEPACAM2_*
do
cd $i/ncbi_dataset/data
/mnt/disk4/BUDDHA/SAMD9-9L/dataformat tsv gene-product --inputfile product_report.jsonl > product_report.dataformat.tsv
cut -f9,18 product_report.dataformat.tsv|tail -n+2|sort -u |sed 's/ /_/g;s/\t/-/g' > scaffold_chr.list
cd /mnt/disk4/BUDDHA/SAMD9-9L
done

mkdir CDK6-HEPACAM2
for order in Pholidota Carnivora Perissodactyla Artiodactyla
do
comm -12 HEPACAM2_"$order"/ncbi_dataset/data/scaffold_chr.list CDK6_"$order"/ncbi_dataset/data/scaffold_chr.list > CDK6-HEPACAM2/"$order".common_scaffold_chr.list
echo $order
done

cd CDK6-HEPACAM2
for order in Pholidota Carnivora Perissodactyla Artiodactyla
do
for i in `cat "$order".common_scaffold_chr.list`
do
sp=`echo $i|cut -f1 -d'-'`
id=`echo $i|cut -f2 -d'-'`
efetch -db nucleotide -format fasta -id $id  > "$sp"_"$id".fasta
done
done

cp ../Homo_sapiens_chr7.fa .
sed -i 's/ .*//;s/\./_/' Homo_sapiens_chr7.fa
conda activate py38
for i in `ls -1 *.fasta|head -1`
do
sed -i 's/ .*//;s/\./_/' $i
query_name=`echo $i|sed 's/\.fasta//g'`
/mnt/disk4/BUDDHA/tools/make_lastz_chains/make_chains.py Homo_sapiens $query_name Homo_sapiens_chr7.fa $i --executor_queuesize 16 
rm -fr TEMP_run.lastz/ TEMP_run.cat/ TEMP_run.fillChain/ TEMP_pslParts/ TEMP_axtChain/ TEMP_psl/ cleanUp.csh DEF make_chains.log make_chains_py_params.json master_script.sh
done

#conda activate py311
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export NXF_VER=22.10.0
chain=$1
query=$2
/mnt/disk4/BUDDHA/tools/TOGA_new/TOGA/toga.py $chain GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.bed12 Gallus_gallus.2bit "$query".2bit --pn TOGA-Gallus_gallus_"$query" --nc /mnt/disk4/BUDDHA/tools/TOGA_new/TOGA/nextflow_config_files/ --cesar_jobs_num 500 --cesar_buckets 3,5,10,15,20 --ces --chain_jobs_num 60 --cesar_mem_limit 150

### TOGA
echo -e "Species_name\tSAMD9L_gene_status\tSAMD9_gene_status" > SAMD9-9L_gene_status.tsv
for i in `ls -1 -d chain-*`
do
cd $i
query_name=`echo $i|cut -f2 -d'-'`
species_name=`echo $query_name|cut -f1,2 -d'_'`
#SAMD9L
/mnt/disk4/BUDDHA/tools/TOGA_new/TOGA/toga.py Homo_sapiens."$query_name".allfilled.chain.gz /mnt/disk4/BUDDHA/SAMD9-9L/Homo_sapiens_SAMD9L.bed12 Homo_sapiens.2bit "$query_name".2bit --pn TOGA_SAMD9L --nc /mnt/disk4/BUDDHA/tools/TOGA_new/TOGA/nextflow_config_files/ --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
id=`tail -n+2 TOGA_SAMD9L/orthology_scores.tsv |sort -nrk3|cut -f1,2|head -1|sed 's/\t/\./g'`
/mnt/disk4/BUDDHA/TOGA-1.1.14/supply/plot_mutations.py --publication_mode_heni /mnt/disk4/BUDDHA/SAMD9-9L/Homo_sapiens_SAMD9L.bed12 TOGA_SAMD9L/inact_mut_data.txt NM_152703.5 "$species_name"-SAMD9L.svg
grep -A1 -w "$id" TOGA_SAMD9L/codon.fasta |tail -2 |sed "s/>.*/>$species_name/g" |sed 's/ //g;s/-//g;s/XXX/NNN/g'> "$species_name"-SAMD9L.fa
cat "$species_name"-SAMD9L.fa >> ../All_species.SAMD9L.fa
SAMD9L_gene_status=`grep -w "$id" TOGA_SAMD9L/loss_summ_data.tsv |cut -f3`
#SAMD9
/mnt/disk4/BUDDHA/tools/TOGA_new/TOGA/toga.py Homo_sapiens."$query_name".allfilled.chain.gz /mnt/disk4/BUDDHA/SAMD9-9L/Homo_sapiens_SAMD9.bed12 Homo_sapiens.2bit "$query_name".2bit --pn TOGA_SAMD9 --nc /mnt/disk4/BUDDHA/tools/TOGA_new/TOGA/nextflow_config_files/ --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
id=`tail -n+2 TOGA_SAMD9/orthology_scores.tsv |sort -nrk3|cut -f1,2|head -1|sed 's/\t/\./g'`
/mnt/disk4/BUDDHA/TOGA-1.1.14/supply/plot_mutations.py --publication_mode_heni /mnt/disk4/BUDDHA/SAMD9-9L/Homo_sapiens_SAMD9.bed12 TOGA_SAMD9/inact_mut_data.txt NM_017654.4 "$species_name"-SAMD9.svg
grep -A1 -w "$id" TOGA_SAMD9/codon.fasta |tail -2 |sed "s/>.*/>$species_name/g" |sed 's/ //g;s/-//g;s/XXX/NNN/g'> "$species_name"-SAMD9.fa
cat "$species_name"-SAMD9.fa >> ../All_species.SAMD9.fa
SAMD9_gene_status=`grep -w "$id" TOGA_SAMD9/loss_summ_data.tsv |cut -f3`
echo -e "$species_name\t$SAMD9L_gene_status\t$SAMD9_gene_status" >> ../SAMD9-9L_gene_status.tsv
cd ..
done

