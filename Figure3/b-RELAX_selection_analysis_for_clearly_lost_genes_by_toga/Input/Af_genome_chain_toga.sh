genome=$1
i=$2
chr_id=$3
species_name=`echo $genome|cut -f1 -d'-'`
blastn -task blastn -evalue 0.05 -db $genome -out $i/blastn_"$i"_cds_"$species_name".bls -num_threads 8 -outfmt 7 -query $i/"$i".fa
/media/morpheus/sagar/BUDDHA/Tools/blast2bed/blast2bed $i/blastn_"$i"_cds_"$species_name".bls
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
/media/morpheus/sagar/BUDDHA/TOGA_new/make_lastz_chains/make_chains.py Antechinus_flavipes $species_name Antechinus_flavipes-"$chr_id".fa "$species_name".30kbslop.fasta --chaining_memory 40 --project_dir chain_Antechinus_flavipes_"$species_name" --kt 
cd chain_Antechinus_flavipes_"$species_name"
mv query.2bit ../"$species_name".2bit
mv target.2bit ../Antechinus_flavipes.2bit
mv Antechinus_flavipes."$species_name".final.chain.gz ../
cd ..
rm -r chain_Antechinus_flavipes_"$species_name"

for chain in Antechinus_flavipes."$species_name".final.chain.gz
do
bed_input="$i"-GCF_016432865.1_AdamAnt_v2_genomic.bed12
id=`cat $bed_input|cut -f4`
/media/morpheus/sagar/BUDDHA/TOGA/toga.py $chain $bed_input Antechinus_flavipes.2bit "$species_name".2bit --pn TOGA_Antechinus_flavipes_"$species_name" --nc /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/ --cesar_bigmem_config /media/morpheus/sagar/BUDDHA/TOGA/nextflow_config_files/cesar_bigmem_config.nf --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
/media/morpheus/sagar/BUDDHA/TOGA/supply/plot_mutations.py --publication_mode_heni $bed_input TOGA_Antechinus_flavipes_"$species_name"/inact_mut_data.txt $id sorted/mutation_plot/"$species_name".svg
done

for toga in TOGA_Antechinus_flavipes_"$species_name"
do
cat $toga/inact_mut_data.txt| sed "s/$id/$species_name/g" >> sorted/"$i".inact_mut_data.txt
orthology_classification=`cat $toga/orthology_classification.tsv|tail -1|cut -f5`; loss_summ=`cat $toga/loss_summ_data.tsv|tail -1|cut -f3`
echo -e "$species_name\t$orthology_classification\t$loss_summ" >> sorted/"$i".orthology.loss_summ.tsv
cat $toga/codon.fasta |tail -2 |sed 's/ //g;s/X/N/g'|sed "s/>.*/>$species_name/g" >> sorted/"$i".codon.fasta
cat $toga/nucleotide.fasta |tail -2 |sed 's/ //g'|sed "s/>.*/>$species_name/g" >> sorted/"$i".nucleotide.fasta
done

cd ..
