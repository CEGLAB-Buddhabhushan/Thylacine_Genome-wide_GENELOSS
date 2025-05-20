# save genome info in genome_info.tsv from table S3
cut -f1,3  genome_info.tsv|sed 's/\t/-/g'|tail -n+2|sed 's/ /_/g' > chr_marsupial.lst
# download genomes ./1_download_data.sh
#create the file chr_marsupial.lst which have acc id and species name seprated by "-" 
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod 777 datasets
for i in `cat chr_marsupial.lst|head -1`
do
acc=`echo $i|cut -f1 -d'-'`
sp=`echo $i|cut -f2 -d'-'`
./datasets download genome accession $acc --filename "$sp"_dataset.zip --include genome --api-key c181c44f2dca87803a9c0f9a24a32f67a608
unzip "$sp"_dataset.zip
genome=`ls -1 ncbi_dataset/data/$acc/*.fna |cut -f4 -d'/'`
mv ncbi_dataset/data/$acc/*.fna "$sp"-"$genome"
rm -r md5sum.txt "$sp"_dataset.zip ncbi_dataset README.md
done

# make balst database and get chr sizes ./2_make_database.sh
ls -1 *.fna > list_of_genomes
for genome in `cat list_of_genomes`
do
makeblastdb -in $genome -out $genome -dbtype nucl
faidx $genome -i chromsizes > "$genome".sizes.genome
done


## get the gtf of tasmanian devil
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/635/505/GCF_902635505.1_mSarHar1.11/GCF_902635505.1_mSarHar1.11_genomic.gtf.gz
gunzip GCF_902635505.1_mSarHar1.11_genomic.gtf.gz
/media/morpheus/sagar/BUDDHA/Tools/ea-utils/clipper/gtf2bed GCF_902635505.1_mSarHar1.11_genomic.gtf > GCF_902635505.1_mSarHar1.11_genomic.bed12


# extact seq , make chain, & run TOGA ./3_21species_parallel_extract_seq_chain_toga.sh
#!/bin/bash

##create transcript.lst file with all gene transcripts needed to be checked
# Initialize variables and environment
i=$1
ulimit -n 16384
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export NXF_VER=22.10.0

# Ensure directories exist and retrieve necessary data
[ ! -d "$i" ] && mkdir -p  $i/sorted/mutation_plot
efetch -db nuccore -id $i -format fasta_cds_na | cut -f1,2 -d'_' | sed 's/>lcl|/>/g' > $i/"$i".fa
length=$(faidx  $i/"$i".fa -i chromsizes | cut -f2)
GC=$(seqkit fx2tab --name --gc  $i/"$i".fa | cut -f2)
echo -e "$i\t$length\t$GC" > $i/"$i".tsv
chr_id=$(grep -w "$i" GCF_902635505.1_mSarHar1.11_genomic.bed12 | cut -f1 | sed 's/_/\./2')
grep -w "$i" GCF_902635505.1_mSarHar1.11_genomic.bed12 > $i/"$i"-GCF_902635505.1_mSarHar1.11_genomic.bed12
echo $chr_id > $i/chr.lst
seqtk subseq Sarcophilus_harrisii-GCF_902635505.1_mSarHar1.11_genomic.fna $i/chr.lst -l 60 > $i/Sarcophilus_harrisii-"$chr_id".fa
sed -i 's/ .*//;s/\./_/' $i/Sarcophilus_harrisii-"$chr_id".fa

# Run parallel jobs for multiple species
./genome_chain_toga.sh Antechinus_flavipes-GCA_016432865.2_AdamAnt_v2_genomic.fna $i $chr_id &
./genome_chain_toga.sh Dasyurus_viverrinus-GCA_020854095.1_UniMelb_DasViv_v1.0_genomic.fna $i $chr_id &
./genome_chain_toga.sh Dromiciops_gliroides-GCA_019393635.1_mDroGli1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Gracilinanus_agilis-GCA_016433145.1_AgileGrace_genomic.fna $i $chr_id &
./genome_chain_toga.sh Lagorchestes_hirsutus-GCA_028533205.1_Lagorchestes_hirsutus_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Macropus_fuliginosus-GCA_028583105.1_mf-2k_genomic.fna $i $chr_id &
./genome_chain_toga.sh Macropus_giganteus-GCA_028627215.1_mg-2k_genomic.fna $i $chr_id &
./genome_chain_toga.sh Macrotis_lagotis-GCA_037893015.1_bilby.v1.9.chrom.fasta_genomic.fna $i $chr_id &
./genome_chain_toga.sh Monodelphis_domestica-GCA_027887165.1_mMonDom1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Myrmecobius_fasciatus-GCA_023553655.1_mMyrfas1.20211206_genomic.fna $i $chr_id &
wait
./genome_chain_toga.sh Notamacropus_eugenii-GCA_028372415.2_mMacEug1.pri_v2_genomic.fna $i $chr_id &
./genome_chain_toga.sh Phalanger_gymnotis-GCA_028646595.1_pg-2k_genomic.fna $i $chr_id &
./genome_chain_toga.sh Potorous_gilbertii-GCA_028658325.1_Potorous_gilbertii_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Pseudocheirus_occidentalis-GCA_028646575.1_Pseudocheirus_occidentalis_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Pseudochirops_corinnae-GCA_028646515.1_Pseudochirops_corinnae_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Pseudochirops_cupreus-GCA_028627135.1_Pseudochirops_cupreus_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Sminthopsis_crassicaudata-GCA_048593235.1_ASM4859323v1_genomic.fna $i $chr_id &
./genome_chain_toga.sh Thylacinus_cynocephalus-GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly_genomic.fna $i $chr_id &
./genome_chain_toga.sh Trichosurus_vulpecula-GCA_011100635.1_mTriVul1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Vombatus_ursinus-GCA_028626985.1_vu-2k_genomic.fna $i $chr_id &
wait
# Continue with the next steps
echo "All jobs for $i completed."



## genome_chain_toga.sh
genome=$1
i=$2
chr_id=$3
species_name=`echo $genome|cut -f1 -d'-'`
blastn -task blastn -evalue 0.05 -db $genome -out $i/blastn_"$i"_cds_"$species_name".6out -num_threads 8 -outfmt '6 sseqid sstart send qseqid evalue sframe qlen qcovs' -query $i/"$i".fa
awk '$8>=75 {print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' $i/blastn_"$i"_cds_"$species_name".6out > $i/blastn_"$i"_cds_"$species_name".bed
awk '$8 >= 75 { strand = ($6 == "-1" ? "-" : ($6 == "1" || $6 == "+1" ? "+" : $6)); print $1, $2, $3, $4, $5, strand }' OFS='\t'

sed -i "/#.*/d;s/$i/$species_name/g" $i/blastn_"$i"_cds_"$species_name".bed
bedtools merge -s -header -d 50000 -i  $i/blastn_"$i"_cds_"$species_name".bed > $i/blastn_"$i"_cds_"$species_name".merge.bed
#faidx $genome -i chromsizes > "$genome".sizes.genome
bedtools slop -i $i/blastn_"$i"_cds_"$species_name".merge.bed -g "$genome".sizes.genome -b 100000 -header > $i/blastn_"$i"_cds_"$species_name".100kbslop.bed
bedtools getfasta -fi $genome -bed $i/blastn_"$i"_cds_"$species_name".100kbslop.bed  > $i/"$species_name".100kbslop.fasta
seqkit rmdup -s  $i/"$species_name".100kbslop.fasta > $i/"$species_name".100kbslop.fas
mv $i/"$species_name".100kbslop.fas $i/"$species_name".100kbslop.fasta
echo $species_name $genome

cd $i

sed -i 's/ .*//;s/\./_/;s/:/_/g;s/-//g;s/[()]//g' "$species_name".100kbslop.fasta
/media/morpheus/sagar/BUDDHA/TOGA_new/make_lastz_chains/make_chains.py Sarcophilus_harrisii $species_name Sarcophilus_harrisii-"$chr_id".fa "$species_name".100kbslop.fasta --chaining_memory 40 --project_dir chain_Sarcophilus_harrisii_"$species_name" --kt 
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
/media/morpheus/sagar/BUDDHA/TOGA/toga.py $chain $bed_input Sarcophilus_harrisii.2bit "$species_name".2bit --pn TOGA_Sarcophilus_harrisii_"$species_name" --nc /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/ --cesar_bigmem_config /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/cesar_bigmem_config.nf --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
/media/morpheus/sagar/BUDDHA/TOGA/supply/plot_mutations.py --publication_mode_heni $bed_input TOGA_Sarcophilus_harrisii_"$species_name"/inact_mut_data.txt $id sorted/mutation_plot/"$species_name".svg
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



################################################################################################################################
genome=$1
i=$2
chr_id=$3
species_name=`echo $genome|cut -f1 -d'-'`

blastn -task blastn -evalue 0.05 -db $genome -out $i/blastn_"$i"_cds_"$species_name".6out -num_threads 8 -outfmt '6 sseqid sstart send qseqid evalue sframe qlen qcovs' -query $i/"$i".fa
awk '$8 >= 75 {
  strand = ($5 == "-1" ? "-" : ($5 == "1" ||$5 == "+1" ? "+" : $5));
  print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"strand
}' $i/blastn_"$i"_cds_"$species_name".6out > $i/blastn_"$i"_cds_"$species_name".bed

awk '$8 >= 75 {
  strand = ($6 == "-1" ? "-" : ($6 == "+1" || $6 == "1" ? "+" : $6));
  start = ($2 < $3 ? $2 : $3);
  end = ($2 > $3 ? $2 : $3);
  print $1"\t"start"\t"end"\t"$4"\t"strand"\t"$5
}' 
sort-bed  $i/blastn_"$i"_cds_"$species_name".bed > $i/blastn_"$i"_cds_"$species_name".sort.bed
mv $i/blastn_"$i"_cds_"$species_name".sort.bed $i/blastn_"$i"_cds_"$species_name".bed
sed -i "/#.*/d;s/$i/$species_name/g" $i/blastn_"$i"_cds_"$species_name".bed
bedtools merge -s -header -d 50000 -i  $i/blastn_"$i"_cds_"$species_name".bed > $i/blastn_"$i"_cds_"$species_name".merge.bed
#faidx $genome -i chromsizes > "$genome".sizes.genome
bedtools slop -i $i/blastn_"$i"_cds_"$species_name".merge.bed -g "$genome".sizes.genome -b 100000 -header > $i/blastn_"$i"_cds_"$species_name".100kbslop.bed
bedtools getfasta -fi $genome -bed $i/blastn_"$i"_cds_"$species_name".100kbslop.bed  > $i/"$species_name".100kbslop.fasta
seqkit rmdup -s  $i/"$species_name".100kbslop.fasta > $i/"$species_name".100kbslop.fas
mv $i/"$species_name".100kbslop.fas $i/"$species_name".100kbslop.fasta
echo $species_name $genome

cd $i

sed -i 's/ .*//;s/\./_/;s/:/_/g;s/-//g;s/[()]//g' "$species_name".100kbslop.fasta
/media/morpheus/sagar/BUDDHA/TOGA_new/make_lastz_chains/make_chains.py Sarcophilus_harrisii $species_name Sarcophilus_harrisii-"$chr_id".fa "$species_name".100kbslop.fasta --chaining_memory 40 --project_dir chain_Sarcophilus_harrisii_"$species_name" --kt 
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
/media/morpheus/sagar/BUDDHA/TOGA/toga.py $chain $bed_input Sarcophilus_harrisii.2bit "$species_name".2bit --pn TOGA_Sarcophilus_harrisii_"$species_name" --nc /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/ --cesar_bigmem_config /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/cesar_bigmem_config.nf --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
/media/morpheus/sagar/BUDDHA/TOGA/supply/plot_mutations.py --publication_mode_heni $bed_input TOGA_Sarcophilus_harrisii_"$species_name"/inact_mut_data.txt $id sorted/mutation_plot/"$species_name".svg
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

