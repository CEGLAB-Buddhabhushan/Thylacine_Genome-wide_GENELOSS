##Extract required columns
mkdir kallisto_output
for i in *.kallisto_out
do
name=`echo $i|cut -f1 -d'.'`
cp $i/abundance.tsv kallisto_output/"$name".abundance.tsv
echo $name
done
cd kallisto_output
cat `ls -1|head -1`|cut -f1 -d$'\t' > counts.tsv
for i in *.abundance.tsv
do
name=`echo $i|cut -f1 -d'.'`
cut -f5 -d$'\t' $i | tail -n +2 > "$name"_headless.tsv
echo -e "$name" | cat - "$name"_headless.tsv > "$name".tsv
rm "$name"_headless.tsv
paste counts.tsv "${name}.tsv" > temp_counts.tsv
mv temp_counts.tsv counts.tsv
done


