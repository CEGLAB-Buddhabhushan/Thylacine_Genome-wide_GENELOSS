#!/bin/bash
echo "TranscriptID,No_Transcripts,Low_quality_protein,Id,Name,Type,Description,Chromosome,ExonCount,ORF status,GC Content,GC Stretch,Missing exon,Deleted exon,Splice site change,START missing,TOGA status,Events supported in SRA,CL ratio"|sed 's/,/\t/g' > UL.info.tsv
for transcriptID in `cat UL.transcript_checked.lst`
do

#transcript_info_5c=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element Id,Name,Description,Chromosome,ExonCount`
Id=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element Id`
Name=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element Name`
Type=`efetch -db nuccore -id $transcriptID -format gb |grep "/mol_type="|cut -f2 -d'='|sed 's/"//g'`
Description=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element Description`
Chromosome=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element Chromosome`
ExonCount=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element ExonCount`

no_transcripts=`grep "$Name" /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/GCF_902635505.1_mSarHar1.11_genomic.gtf|awk '$3=="transcript"'|wc -l`
orf_status=`perl /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/ORFvalid.pl $transcriptID/"$transcriptID".fa`
if echo "$orf_status" | grep -q "1 out of 1 sequences validated as ORFs."; then
    orf_val="yes"
else
    orf_val="no"
fi
GC_Content=`seqkit fx2tab --name --gc $transcriptID/"$transcriptID".fa | awk '{print $2}'`
GC_Stretch=`perl /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/GC_Stretch_finder.pl $transcriptID/"$transcriptID".fa |awk '{print $5}'`

# Check for exon missing and deleted
if grep -q "Missing exon" $transcriptID/"$transcriptID".inact_mut_data.txt ; then
exon_missing="yes"
else
exon_missing="no"
fi

if grep -q "Deleted exon" $transcriptID/"$transcriptID".inact_mut_data.txt ; then
exon_deleted="yes"
else
exon_deleted="no"
fi

if grep -q "SSM" $transcriptID/"$transcriptID".inact_mut_data.txt ; then
splice_site_change="yes"
else
splice_site_change="no"
fi

if grep -q "START_MISSING" $transcriptID/"$transcriptID".inact_mut_data.txt ; then
start_missing="yes"
else
start_missing="no"
fi
if efetch -db nuccore -id $transcriptID -format gb|grep -q "LOW QUALITY PROTEIN"; then
low_quality_protein="yes"
else
low_quality_protein="no"
fi
sra_confirm=`awk 'BEGIN{total=0; yes_count=0} {total++; if($5 == "yes") yes_count++} END {print yes_count "/" total, yes_count/total}' $transcriptID/"$transcriptID".codon_positions.events.TOGA_confirm.tsv|sed 's/ /\t/g'`
toga_status=`grep -w "$transcriptID" /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA.loss_summ_data.tsv|cut -f2`
echo -e "$transcriptID\t$no_transcripts\t$low_quality_protein\t$Id\t$Name\t$Type\t$Description\t$Chromosome\t$ExonCount\t$orf_val\t$GC_Content\t$GC_Stretch\t$exon_missing\t$exon_deleted\t$splice_site_change\t$start_missing\t$toga_status\t$sra_confirm" >> UL.info.tsv
echo $transcriptID
done

