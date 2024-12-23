#create the file chr_marsupial.lst which have acc id and species name seprated by "-" 
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod 777 datasets
for i in `cat chr_marsupial.lst`
do
acc=`echo $i|cut -f1 -d'-'`
sp=`echo $i|cut -f2 -d'-'`
./datasets download genome accession $acc --filename "$sp"_dataset.zip --include genome
unzip "$sp"_dataset.zip
genome=`ls -1 ncbi_dataset/data/$acc/*.fna |cut -f4 -d'/'`
mv ncbi_dataset/data/$acc/*.fna "$sp"-"$genome"
rm -r md5sum.txt "$sp"_dataset.zip ncbi_dataset README.md
done

