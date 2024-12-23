#!/bin/bash
### provide input like ./Confirm_Loss.sh transcriptID
transcriptID=$1

# Download CDS
efetch -db nuccore -id $transcriptID -format fasta_cds_na|cut -f1,2 -d'_'|sed 's/>lcl|/>/g' > "$transcriptID".fa
# Transcript info and preliminary assesment

transcript_info_8c=`esearch -db gene -query $transcriptID | esummary | xtract -pattern DocumentSummary -element Id,Name,Description,Chromosome,ChrAccVer,ChrStart,ChrStop,ExonCount`
orf_status=`perl /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/ORFvalid.pl "$transcriptID".fa`
if echo "$orf_status" | grep -q "1 out of 1 sequences validated as ORFs."; then
    orf_val="yes"
else
    orf_val="no"
fi
GC_Content=`seqkit fx2tab --name --gc "$transcriptID".fa | awk '{print $2}'`
GC_Stretch=`perl /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/GC_Stretch_finder.pl "$transcriptID".fa |awk '{print $5}'`

# Convert TOGA evenets
python /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/make_codon_position_bed.py $transcriptID
grep -w "$transcriptID" /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA.inact_mut_data.txt > "$transcriptID".inact_mut_data.txt
# Check for exon missing and deleted
if grep -q "Missing exon" "$transcriptID".inact_mut_data.txt ; then
exon_missing="yes"
else
exon_missing="no"
fi

if grep -q "Deleted exon" "$transcriptID".inact_mut_data.txt ; then
exon_deleted="yes"
else
exon_deleted="no"
fi

if grep -q "SSM" "$transcriptID".inact_mut_data.txt ; then
splice_site_change="yes"
else
splice_site_change="no"
fi

if grep -q "START_MISSING" "$transcriptID".inact_mut_data.txt ; then
start_missing="yes"
else
start_missing="no"
fi

cat "$transcriptID".inact_mut_data.txt |cut -f3-6|grep -v "GENE:"|awk 'NF == 4 && $1!=0 {print "Exon"$1":codon"$2":"$3":"$4}' > "$transcriptID".events.bed

bash /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/generate_codon_bed_for_cds.sh $transcriptID
for i in `cat "$transcriptID".events.bed`; do pos=`echo $i|cut -f2 -d':'`; grep -w "$pos" "$transcriptID".cds_codons.bed|sed "s/$pos/$i/g"|awk '{print $0"\t"$2"\t"$3"\t"255","0","0}'  >> "$transcriptID".codon_positions.events.bed; done

# BLASTn in DNAseqDB
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db /media/morpheus/disk1/BUDDHA_merged_data/ThyCyn_db/SRR5055303-6.blastDB.fa -out "$transcriptID".blastn.DNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query "$transcriptID".fa
sed -i "s/Query_1/$transcriptID/g" "$transcriptID".blastn.DNAseqDB.sam
samtools view -bhS "$transcriptID".blastn.DNAseqDB.sam > "$transcriptID".blastn.DNAseqDB.bam
samtools sort "$transcriptID".blastn.DNAseqDB.bam -o "$transcriptID".blastn.DNAseqDB.sorted.bam
samtools index "$transcriptID".blastn.DNAseqDB.sorted.bam

# BLASTn in RNAseqDB
blastn -task blastn -evalue 0.01 -max_target_seqs 5000 -db /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/mapping/miRNA_seq/miRNAblastDB/SRR23147611-6.miRNAblastDB.fa -out "$transcriptID".blastn.miRNAseqDB.sam -num_threads 32 -outfmt '17 SQ'  -query "$transcriptID".fa
sed -i "s/Query_1/$transcriptID/g" "$transcriptID".blastn.miRNAseqDB.sam
samtools view -bhS "$transcriptID".blastn.miRNAseqDB.sam > "$transcriptID".blastn.miRNAseqDB.bam
samtools sort "$transcriptID".blastn.miRNAseqDB.bam -o "$transcriptID".blastn.miRNAseqDB.sorted.bam
samtools index "$transcriptID".blastn.miRNAseqDB.sorted.bam

# Pilon correction
java -Xmx300g -jar /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/pilon-1.24.jar --genome "$transcriptID".fa --frags "$transcriptID".blastn.DNAseqDB.sorted.bam --output "$transcriptID".pilon.DNAseq --fix all --changes

# create IGV report
create_report "$transcriptID".codon_positions.events.bed --standalone --fasta "$transcriptID".fa --tracks "$transcriptID".blastn.DNAseqDB.sorted.bam "$transcriptID".codon_positions.events.bed --output "$transcriptID".BLASTn_IGV.html --info-columns Chromosome Start_position End_position Event score strand thickStart thickEnd itemRgb --translate-sequence-track

# Short-read verification using bam-readcount
awk '{print $1"\t"$2"\t"$3"\t"$4; if ($4 ~ /->/) {split($4, a, ":"); split(a[4], change, "->"); start = $2 + 1; for (i = 1; i <= length(change[1]); i++) {ref_base = substr(change[1], i, 1); alt_base = substr(change[2], i, 1); print $1"\t"start"\t"ref_base"\t"alt_base; start++}}}' "$transcriptID".codon_positions.events.bed | grep -v "STOP" | awk '$3!=$4' | awk '{if ($4 ~ /Exon/) {split($4, a, ":"); print $1"\t"$2"-"$3"\t"a[3]"\t"a[4]} else {if ($2 == $3) {print $1"\t"$2"-"$2"\t"$3"\t"$4} else {print $1"\t"$2"-"$2"\t"$3"\t"$4}}}'|sed 's/\t/,/g' > "$transcriptID".codon_positions.events.brc.csv
for i in `cat "$transcriptID".codon_positions.events.brc.csv`; do
  
  file_name=`echo $i | sed 's/,/_/g'`
  chr=`echo $i | cut -f1 -d','`
  start_end=`echo $i | cut -f2 -d','`
  mutation_type=`echo $i | cut -f3 -d','`   # Get the 3rd column (mutation type)
  ref_base=`echo $i | cut -f4 -d','`        # Get the 4th column (ref base)
  query_base=`echo $i | cut -f5 -d','`      # Get the 5th column (query base)

  if [[ "$mutation_type" != "FS_INS" && "$mutation_type" != "FS_DEL" ]]; then
    # Normal processing if not FS_INS or FS_DEL
    bam-readcount -w0 -f "$transcriptID".fa "$transcriptID".blastn.DNAseqDB.sorted.bam $chr:$start_end > "$transcriptID"."$chr":"$start_end".brc.tsv
    python /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/parse_brc.py "$transcriptID"."$chr":"$start_end".brc.tsv > "$transcriptID"."$chr":"$start_end".parse_brc.tsv
    event_status=$(awk -v chr="$chr" -v ref_base="$ref_base" -v query_base="$query_base" 'BEGIN {found=0} $7 >= 5 && !found {print "yes"; found=1}' "$transcriptID"."$chr":"$start_end".parse_brc.tsv)
    echo -e "$i,$event_status"

  else
    # Check for + sign in 4th column if FS_INS or - sign if FS_DEL
    bam-readcount -w0 -f "$transcriptID".fa "$transcriptID".blastn.DNAseqDB.sorted.bam $chr:$start_end > "$transcriptID"."$chr":"$start_end".brc.tsv
    python /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/parse_brc.py "$transcriptID"."$chr":"$start_end".brc.tsv > "$transcriptID"."$chr":"$start_end".parse_brc.tsv
    
    if [[ "$mutation_type" == "FS_INS" ]]; then
      # Check for "+" in 4th column of parse_brc.tsv
      if cut -f4 "$transcriptID"."$chr":"$start_end".parse_brc.tsv | grep -q "^+" ; then
        echo -e "$i,yes"
      else
        echo -e "$i,No"
      fi

    elif [[ "$mutation_type" == "FS_DEL" ]]; then
      # Check for "-" in 4th column of parse_brc.tsv
      if cut -f4 "$transcriptID"."$chr":"$start_end".parse_brc.tsv | grep -q "^-"; then
        echo -e "$i,yes"
      else
        echo -e "$i,No"
      fi
    fi
  fi
done|sed 's/,/\t/g'> "$transcriptID".codon_positions.events.TOGA_confirm.tsv
sra_confirm=`awk 'BEGIN{total=0; yes_count=0} {total++; if($5 == "yes") yes_count++} END {print yes_count "/" total, yes_count/total}' "$transcriptID".codon_positions.events.TOGA_confirm.tsv|sed 's/ /\t/g'`
toga_status=`grep -w "$transcriptID" /media/morpheus/sagar/BUDDHA/Tasmanian_wolf/Chr_wise/Chromosomes/Final_verification/TOGA.loss_summ_data.tsv|cut -f2`
echo "TranscriptID,Id,Name,Description,Chromosome,ChrAccVer,ChrStart,ChrStop,ExonCount,ORF status,GC Content,GC Stretch,Missing exon,Deleted exon,Splice site change,START missing,TOGA status,Events supported in SRA,CL ratio"|sed 's/,/\t/g' > "$transcriptID".info.tsv
echo -e "$transcriptID\t$transcript_info_8c\t$orf_val\t$GC_Content\t$GC_Stretch\t$exon_missing\t$exon_deleted\t$splice_site_change\t$start_missing\t$toga_status\t$sra_confirm" >> "$transcriptID".info.tsv

mkdir $transcriptID
mv "$transcriptID".*  $transcriptID
