## Chain files and TOGA output are uploaded to "https://github.com/Buddha358/Thylacine_Genome-wide_GENELOSS-Chromosome_wise_Chains_and_TOGA_results.git" due to large file sizes. ##
Further, chain files split into multiple chain files using the following script:
```
for i in *.Tasmanian_wolf.final.chain.gz
do
gunzip $i
name=`echo $i|sed 's/\.gz//g'`
split -b 5000M $name split_chain."$name".
rm $name
for j in split_chain."$name".*
do
gzip $j
done
echo $i
done
```
