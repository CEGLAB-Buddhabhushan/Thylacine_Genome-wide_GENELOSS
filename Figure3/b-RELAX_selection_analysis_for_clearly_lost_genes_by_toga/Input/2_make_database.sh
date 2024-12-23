ls -1 *.fna > list_of_genomes
for genome in `cat list_of_genomes`
do
makeblastdb -in $genome -out $genome -dbtype nucl
faidx $genome -i chromsizes > "$genome".sizes.genome
done

