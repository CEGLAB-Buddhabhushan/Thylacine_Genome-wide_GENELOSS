for i in `cat Trehalase.Species.lst`
do
sp=`echo $i|sed 's/_/ /g'`
datasets download genome taxon "$sp" --include genome --api-key c181c44f2dca87803a9c0f9a24a32f67a608 --reference --filename "$i".zip --no-progressbar
echo $sp
done


for i in `ls -1 *.zip`
do
sp=`echo $i |sed 's/\.zip//g'`
unzip $i
cd ncbi_dataset/data/*/
genome=`ls -1 *.fna|head -1`
makeblastdb -in $genome -out $genome -dbtype nucl
blastn -task blastn  -evalue 0.05 -db $genome -out blastn.Homo_sapiens_CDK6."$genome".6out -num_threads 64 -outfmt '6 qseqid sseqid evalue qlen qcovs' -query /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/Homo_sapiens_CDK6.fa
blastn -task blastn  -evalue 0.05 -db $genome -out blastn.Homo_sapiens_SAMD9."$genome".6out -num_threads 64 -outfmt '6 qseqid sseqid evalue qlen qcovs' -query /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/Homo_sapiens_SAMD9.fa
blastn -task blastn  -evalue 0.05 -db $genome -out blastn.Homo_sapiens_SAMD9L."$genome".6out -num_threads 64 -outfmt '6 qseqid sseqid evalue qlen qcovs' -query /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/Homo_sapiens_SAMD9L.fa
cut -f2 blastn.Homo_sapiens_CDK6."$genome".6out |sort -u > CDK6.lst
cut -f2 blastn.Homo_sapiens_SAMD9."$genome".6out |sort -u > SAMD9.lst
cut -f2 blastn.Homo_sapiens_SAMD9L."$genome".6out |sort -u > SAMD9L.lst
comm -12 CDK6.lst SAMD9L.lst|comm -12 - SAMD9.lst|sort -u  > /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/"$sp".CDK6_SAMD9_SAMD9L.lst

seqtk subseq $genome /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/"$sp".CDK6_SAMD9_SAMD9L.lst > /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/"$sp".CDK6_SAMD9_SAMD9L.fa
cd /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/
rm -r ncbi_dataset/ md5sum.txt README.md
echo $sp
done

bash chain.sh 


#!/bin/bash

counter=0

for j in *.CDK6_SAMD9_SAMD9L.fa; do
    (
        # Clean up FASTA headers and filename
        sed -i 's/ .*//;s/\./_/' "$j"
        query_name=$(echo "$j" | awk -F'.' '{print $1}')

        # Run the chain generation script
        /mnt/disk4/BUDDHA/tools/make_lastz_chains/make_chains.py Homo_sapiens "$query_name" Homo_sapiens_chr7.fa "$j" --executor_queuesize 16 --project_dir Chain-"$query_name"

        # Clean up unnecessary temp files
        cd Chain-"$query_name"
        rm -fr TEMP_run.lastz/ TEMP_run.cat/ TEMP_run.fillChain/ TEMP_pslParts/ TEMP_axtChain/ TEMP_psl/ cleanUp.csh DEF make_chains.log make_chains_py_params.json master_script.sh
        cd ..
    ) &

    ((counter++))
    if (( counter % 15 == 0 )); then
        wait  # Wait for 15 jobs to finish
    fi
done

wait  # Wait for any remaining jobs


mkdir -p /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/sorted/mutation_plot/
for i in `ls -1 -d Chain-*`
do
cd $i
query_name=$(echo "$i" | awk -F'-' '{print $2}')
/mnt/disk4/BUDDHA/TOGA-1.1.14/toga.py Homo_sapiens."$query_name".allfilled.chain.gz /mnt/disk4/BUDDHA/SAMD9-9L/Homo_sapiens_SAMD9.bed12 Homo_sapiens.2bit "$query_name".2bit --pn TOGA-"$query_name" --nc /mnt/disk4/BUDDHA/TOGA-1.1.14/nextflow_config_files/ --cesar_bigmem_config /mnt/disk4/BUDDHA/TOGA-1.1.14/nextflow_config_files/cesar_bigmem_config.nf --cesar_jobs_num 500 --cesar_buckets 3,5,25,50 --ces --kt --chain_jobs_num 60
/mnt/disk4/BUDDHA/TOGA-1.1.14/supply/plot_mutations.py --publication_mode_heni /mnt/disk4/BUDDHA/SAMD9-9L/Homo_sapiens_SAMD9.bed12 TOGA-"$query_name"/inact_mut_data.txt NM_017654.4 /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/sorted/mutation_plot/"$query_name".svg

cat TOGA-"$query_name"/inact_mut_data.txt| sed "s/NM_017654.4/$query_name/g" >> /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/sorted/SAMD9.inact_mut_data.txt
orthology_classification=`cat TOGA-"$query_name"/orthology_classification.tsv|tail -1|cut -f5`
loss_summ=`cat TOGA-"$query_name"/loss_summ_data.tsv|tail -1|cut -f3`
echo -e "$query_name\t$orthology_classification\t$loss_summ" >> /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/sorted/SAMD9.orthology.loss_summ.tsv
cat TOGA-"$query_name"/codon.fasta |tail -2 |sed 's/ //g;s/X/N/g'|sed "s/>.*/>$query_name/g" >> /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/sorted/"$clade".codon.fasta
cat TOGA-"$query_name"/nucleotide.fasta |tail -2 |sed 's/ //g'|sed "s/>.*/>$query_name/g" >> /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/sorted/"$clade".nucleotide.fasta
cd /mnt/disk4/BUDDHA/SAMD9-9L/Diet_relation/
done


mv sorted SAMD9_sorted
cd SAMD9_sorted
echo -e "species_name,SAMD9_gene_status" > SAMD9.TOGA_inferred_gene_status.csv
awk '$3=="I" || $3=="L"' SAMD9.orthology.loss_summ.tsv |cut -f1,3|sed 's/\t/,/g' >> SAMD9.TOGA_inferred_gene_status.csv
cut -f1 -d',' SAMD9.TOGA_inferred_gene_status.csv > Species.lst

head -n 1 Trehalase_data.csv | awk '{print "Species_name", "SAMD9_gene_status", $0}'|sed 's/ /,/g' > final_output.csv


for i in `tail -n+2 Species.lst`
do
sp=`grep -w "$i" SAMD9.TOGA_inferred_gene_status.csv|cut -f1 -d','`
SAMD9_status=`grep -w "$i" SAMD9.TOGA_inferred_gene_status.csv|cut -f2 -d','`
info=`grep -w "$i" Trehalase_data.csv`
echo -e "$sp,$SAMD9_status,$info" >> final_output.csv
echo $sp
done

