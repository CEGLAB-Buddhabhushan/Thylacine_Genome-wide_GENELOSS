curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat

for i in SAMD9-54809
do
gene_id=`echo $i|cut -f2 -d'-'`
./datasets download gene gene-id $gene_id --ortholog all --no-progressbar --include gene,rna,cds,product-report,protein --filename "$i".zip
rm -rf $i
unzip -q -o "$i".zip -d $i
rm "$i".zip
echo $i
done


for i in SAMD9-54809
do
cd $i/ncbi_dataset/data
rm *.tsv *.gb cds.seqtksubseq.fa Species_TaxonomicID.lst *.fa
/home/ceglab358/BUDDHA/upload-GITHUB_Thylacine_Genome-wide_GENELOSS/R1_RSPB/SAMD9_ncbi_datasets/dataformat tsv gene-product --inputfile product_report.jsonl > product_report.dataformat.tsv
head -1 product_report.dataformat.tsv|cut -f1,3,7-11,16,17,18,33 > LOW_QUALITY_PROTEIN.tsv
grep "LOW QUALITY PROTEIN" product_report.dataformat.tsv |cut -f1,3,7-11,16,17,18,33|sort -u >> LOW_QUALITY_PROTEIN.tsv
cut -f8 product_report.dataformat.tsv| tail  -n+2 |sort -u > Species_TaxonomicID.lst
echo -e "Acc\tSpecies_name\tInternal_stop_codon\tframeshifts\tAnnotated_introns\tCoverage\tGC_content\tNote_comment" > LOW_QUALITY_PROTEIN.efetch.tsv
for acc in `cut -f6 LOW_QUALITY_PROTEIN.tsv|tail -n+2`
do
echo $acc >>  /home/ceglab358/BUDDHA/upload-GITHUB_Thylacine_Genome-wide_GENELOSS/R1_RSPB/SAMD9_ncbi_datasets/Extract.log.out
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=gb&retmode=txt&api_key=61131f45442b718e60396eb45ce73b65bc08">> "$acc".gb
organism=`grep "/organism" "$acc".gb|sed 's/organism="//g;s/[/]//g;s/"//g'`
Note=`awk '/\/note="/ {flag=1; print; next} flag {if (/"/) {flag=0; print} else {print}}' "$acc".gb | sed 's/"//g;s/note=//g;s/[/]//g' | tr '\n' ' ' | sed 's/  */ /g;s/^ //;s/ $//'`
#check internal Stop
if grep -A3 "##RefSeq-Attributes-START##" "$acc".gb | grep -v "##" | sed 's/^ *//' | grep -q "internal stop codons"; then
    internal_stop="Yes"
else
    internal_stop="No"
fi
#Check for "frameshifts"
if grep -A3 "##RefSeq-Attributes-START##" "$acc".gb | grep -v "##" | sed 's/^ *//' | grep -q "frameshifts"; then
    frameshifts="Yes"
else
    frameshifts="No"
fi
#Check "support for all annotated introns"
if echo $Note|grep -q "support for all annotated introns"; then
    annotated_introns="Yes"
else
    annotated_introns="No"
fi
#Get coverage
coverage=`echo $Note|grep -oP '\d+% coverage'|sed 's/ coverage//g'`
#GC content
echo $acc > "$acc".seqtksubseq.list
cut -f1 -d':' cds.fna > cds.seqtksubseq.fa
seqtk subseq cds.seqtksubseq.fa "$acc".seqtksubseq.list > "$acc".seqtksubseq.fa
GC_content=`seqkit fx2tab --name --gc "$acc".seqtksubseq.fa|cut -f2 `

echo -e "$acc\t$organism\t$internal_stop\t$frameshifts\t$annotated_introns\t$coverage\t$GC_content\t$Note" >> LOW_QUALITY_PROTEIN.efetch.tsv
rm "$acc".seqtksubseq.list "$acc".seqtksubseq.fa
done
rm cds.seqtksubseq.fa 
cd /home/ceglab358/BUDDHA/upload-GITHUB_Thylacine_Genome-wide_GENELOSS/R1_RSPB/SAMD9_ncbi_datasets
echo "$i" >> /home/ceglab358/BUDDHA/upload-GITHUB_Thylacine_Genome-wide_GENELOSS/R1_RSPB/SAMD9_ncbi_datasets/Extract.log.out
echo $i
done
