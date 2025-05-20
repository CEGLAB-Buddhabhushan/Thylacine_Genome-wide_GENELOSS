
head -1 S.harrisii.info.plot.tsv | awk '{print $0 "\tParalogs\tNo_of_paralogs"}' > Paralog_info.S.harrisii.info.plot.tsv

for i in $(cut -f4 Tasmanian_devil_paralogue_mart_export.txt | sort -u); do     nu=$(grep -w "$i" Tasmanian_devil_paralogue_mart_export.txt | wc -l);     awk -v gene="$i" -v no="$nu" '$1 == gene {print $0, "yes", no}' S.harrisii.info.plot.tsv; done > Paralog_info.S.harrisii.info.plot.tsv


head -1 Tasmanian_devil_paralogue_mart_export.txt | awk '{print $0 "\tTOGA_status"}' > TOGA_info.Tasmanian_devil_paralogue_mart_export.tsv

for i in $(cut -f4 Tasmanian_devil_paralogue_mart_export.txt | sort -u); do TOGA_status=`grep -w "$i" S.harrisii.info.plot.tsv|cut -f9|head -1`; awk -v gene="$i" '$4 == gene {print $0, TOGA_status}' Tasmanian_devil_paralogue_mart_export.txt; done > TOGA_info.Tasmanian_devil_paralogue_mart_export.tsv


# First, write the header
head -1 Tasmanian_devil_paralogue_mart_export.txt | awk '{print $0 "\tTOGA_status"}' > TOGA_info.Tasmanian_devil_paralogue_mart_export.tsv

# Then annotate each row
for i in $(cut -f4 Tasmanian_devil_paralogue_mart_export.txt | sort -u); do
    TOGA_status=$(awk -v gene="$i" '$1 == gene {print $9; exit}' S.harrisii.info.plot.tsv)
    awk -v gene="$i" -v status="$TOGA_status" '$4 == gene {print $0 "\t" status}' Tasmanian_devil_paralogue_mart_export.txt
done >> TOGA_info.Tasmanian_devil_paralogue_mart_export.tsv
awk 'NR==1 || $NF != ""' TOGA_info.Tasmanian_devil_paralogue_mart_export.tsv > TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv


# Write header
head -1 Tasmanian_devil_paralogue_mart_export.txt | awk '{print $0 "\tTOGA_status"}' > TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv

# Loop to annotate and filter
for i in $(cut -f4 Tasmanian_devil_paralogue_mart_export.txt | sort -u); do
    TOGA_status=$(awk -v gene="$i" '$1 == gene {print $9; exit}' S.harrisii.info.plot.tsv)
    if [[ -n "$TOGA_status" ]]; then
        awk -v gene="$i" -v status="$TOGA_status" '$4 == gene {print $0 "\t" status}' Tasmanian_devil_paralogue_mart_export.txt
    fi
done >> TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv

cut -f19 TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv | head

