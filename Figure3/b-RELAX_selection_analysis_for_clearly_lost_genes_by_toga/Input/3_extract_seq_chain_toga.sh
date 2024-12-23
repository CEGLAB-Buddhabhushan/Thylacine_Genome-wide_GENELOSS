##create transcript.lst file with all gene transcript needed to be checked
#download Sarcophilus_harrisii gtf nad create bed12 
i=$1
ulimit -n 16384
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export NXF_VER=22.10.0
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/635/505/GCF_902635505.1_mSarHar1.11/GCF_902635505.1_mSarHar1.11_genomic.gtf.gz
#gunzip GCF_902635505.1_mSarHar1.11_genomic.gtf.gz
#/mnt/disk4/BUDDHA/tools/ea-utils/clipper/gtf2bed GCF_902635505.1_mSarHar1.11_genomic.gtf |sed 's/\./_/1' > GCF_902635505.1_mSarHar1.11_genomic.bed12

#for i in `cat transcript.lst`
#do
[ ! -d "$i" ] && mkdir -p  $i/sorted/mutation_plot
efetch -db nuccore -id $i -format fasta_cds_na|cut -f1,2 -d'_'|sed 's/>lcl|/>/g' > $i/"$i".fa
length=`faidx  $i/"$i".fa -i chromsizes |cut -f2`
GC=`seqkit fx2tab --name --gc  $i/"$i".fa|cut -f2`
echo -e "$i\t$length\t$GC" > $i/"$i".tsv
chr_id=`grep -w "$i" GCF_902635505.1_mSarHar1.11_genomic.bed12|cut -f1|sed 's/_/\./2'`
grep -w "$i" GCF_902635505.1_mSarHar1.11_genomic.bed12 > $i/"$i"-GCF_902635505.1_mSarHar1.11_genomic.bed12
echo $chr_id > $i/chr.lst
seqtk subseq Sarcophilus_harrisii-GCF_902635505.1_mSarHar1.11_genomic.fna $i/chr.lst -l 60 > $i/Sarcophilus_harrisii-"$chr_id".fa
sed -i 's/ .*//;s/\./_/' $i/Sarcophilus_harrisii-"$chr_id".fa

for genome in `cat list_of_genomes`
do
species_name=`echo $genome|cut -f1 -d'-'`
blastn -task blastn -evalue 0.05 -db $genome -out $i/blastn_"$i"_cds_"$species_name".bls -num_threads 8 -outfmt 7 -query $i/"$i".fa
blast2bed $i/blastn_"$i"_cds_"$species_name".bls
sort-bed  $i/blastn_"$i"_cds_"$species_name".bed > $i/blastn_"$i"_cds_"$species_name".sort.bed
mv $i/blastn_"$i"_cds_"$species_name".sort.bed $i/blastn_"$i"_cds_"$species_name".bed
sed -i "/#.*/d;s/$i/$species_name/g" $i/blastn_"$i"_cds_"$species_name".bed
bedtools merge -s -header -d 10000 -i  $i/blastn_"$i"_cds_"$species_name".bed > $i/blastn_"$i"_cds_"$species_name".merge.bed
#faidx $genome -i chromsizes > "$genome".sizes.genome
bedtools slop -i $i/blastn_"$i"_cds_"$species_name".merge.bed -g "$genome".sizes.genome -b 30000 -header > $i/blastn_"$i"_cds_"$species_name".30kbslop.bed
bedtools getfasta -fi $genome -bed $i/blastn_"$i"_cds_"$species_name".30kbslop.bed  > $i/"$species_name".30kbslop.fasta
seqkit rmdup -s  $i/"$species_name".30kbslop.fasta > $i/"$species_name".30kbslop.fas
mv $i/"$species_name".30kbslop.fas $i/"$species_name".30kbslop.fasta
echo $species_name $genome

cd $i

sed -i 's/ .*//;s/\./_/;s/:/_/g;s/-//g;s/[()]//g' "$species_name".30kbslop.fasta
/mnt/disk4/BUDDHA/tools/TOGA_new/make_lastz_chains/make_chains.py Sarcophilus_harrisii $species_name Sarcophilus_harrisii-"$chr_id".fa "$species_name".30kbslop.fasta --chaining_memory 40 --project_dir chain_Sarcophilus_harrisii_"$species_name" --kt 
cd chain_Sarcophilus_harrisii_"$species_name"
mv query.2bit ../"$species_name".2bit
mv target.2bit ../Sarcophilus_harrisii.2bit
mv Sarcophilus_harrisii."$species_name".final.chain.gz ../
cd ..
rm -r chain_Sarcophilus_harrisii_"$species_name"

for chain in Sarcophilus_harrisii."$species_name".final.chain.gz
do
bed_input="$i"-GCF_902635505.1_mSarHar1.11_genomic.bed12
id=`cat $bed_input|cut -f4`
/mnt/disk4/BUDDHA/TOGA-1.1.14/toga.py $chain $bed_input Sarcophilus_harrisii.2bit "$species_name".2bit --pn TOGA_Sarcophilus_harrisii_"$species_name" --nc /mnt/disk4/BUDDHA/TOGA-1.1.14/nextflow_config_files/ --cesar_bigmem_config /mnt/disk4/BUDDHA/TOGA-1.1.14/nextflow_config_files/cesar_bigmem_config.nf --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
/mnt/disk4/BUDDHA/TOGA-1.1.14/supply/plot_mutations.py --publication_mode_heni $bed_input TOGA_Sarcophilus_harrisii_"$species_name"/inact_mut_data.txt $id sorted/mutation_plot/"$species_name".svg
done

for toga in TOGA_Sarcophilus_harrisii_"$species_name"
do
cat $toga/inact_mut_data.txt| sed "s/$id/$species_name/g" >> sorted/"$i".inact_mut_data.txt
orthology_classification=`cat $toga/orthology_classification.tsv|tail -1|cut -f5`; loss_summ=`cat $toga/loss_summ_data.tsv|tail -1|cut -f3`
echo -e "$species_name\t$orthology_classification\t$loss_summ" >> sorted/"$i".orthology.loss_summ.tsv
cat $toga/codon.fasta |tail -2 |sed 's/ //g;s/X/N/g'|sed "s/>.*/>$species_name/g" >> sorted/"$i".codon.fasta
cat $toga/nucleotide.fasta |tail -2 |sed 's/ //g'|sed "s/>.*/>$species_name/g" >> sorted/"$i".nucleotide.fasta
done

cd ..

done
