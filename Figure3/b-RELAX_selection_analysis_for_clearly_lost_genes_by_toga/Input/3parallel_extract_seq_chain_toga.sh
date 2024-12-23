#!/bin/bash

##create transcript.lst file with all gene transcripts needed to be checked
# Initialize variables and environment
i=$1
ulimit -n 16384
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export NXF_VER=22.10.0

# Ensure directories exist and retrieve necessary data
[ ! -d "$i" ] && mkdir -p  $i/sorted/mutation_plot
efetch -db nuccore -id $i -format fasta_cds_na | cut -f1,2 -d'_' | sed 's/>lcl|/>/g' > $i/"$i".fa
length=$(faidx  $i/"$i".fa -i chromsizes | cut -f2)
GC=$(seqkit fx2tab --name --gc  $i/"$i".fa | cut -f2)
echo -e "$i\t$length\t$GC" > $i/"$i".tsv
chr_id=$(grep -w "$i" GCF_902635505.1_mSarHar1.11_genomic.bed12 | cut -f1 | sed 's/_/\./2')
grep -w "$i" GCF_902635505.1_mSarHar1.11_genomic.bed12 > $i/"$i"-GCF_902635505.1_mSarHar1.11_genomic.bed12
echo $chr_id > $i/chr.lst
seqtk subseq Sarcophilus_harrisii-GCF_902635505.1_mSarHar1.11_genomic.fna $i/chr.lst -l 60 > $i/Sarcophilus_harrisii-"$chr_id".fa
sed -i 's/ .*//;s/\./_/' $i/Sarcophilus_harrisii-"$chr_id".fa

# Run parallel jobs for multiple species
./genome_chain_toga.sh Antechinus_flavipes-GCA_016432865.2_AdamAnt_v2_genomic.fna $i $chr_id &
./genome_chain_toga.sh Dasyurus_viverrinus-GCA_020854095.1_UniMelb_DasViv_v1.0_genomic.fna $i $chr_id &
./genome_chain_toga.sh Dromiciops_gliroides-GCA_019393635.1_mDroGli1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Gracilinanus_agilis-GCA_016433145.1_AgileGrace_genomic.fna $i $chr_id &
./genome_chain_toga.sh Lagorchestes_hirsutus-GCA_028533205.1_Lagorchestes_hirsutus_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Macropus_fuliginosus-GCA_028583105.1_mf-2k_genomic.fna $i $chr_id &
./genome_chain_toga.sh Macropus_giganteus-GCA_028627215.1_mg-2k_genomic.fna $i $chr_id &
./genome_chain_toga.sh Macrotis_lagotis-GCA_037893015.1_bilby.v1.9.chrom.fasta_genomic.fna $i $chr_id &
./genome_chain_toga.sh Monodelphis_domestica-GCA_027887165.1_mMonDom1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Notamacropus_eugenii-GCA_028372415.1_mMacEug1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Phalanger_gymnotis-GCA_028646595.1_pg-2k_genomic.fna $i $chr_id &
./genome_chain_toga.sh Potorous_gilbertii-GCA_028658325.1_Potorous_gilbertii_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Pseudocheirus_occidentalis-GCA_028646575.1_Pseudocheirus_occidentalis_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Pseudochirops_corinnae-GCA_028646515.1_Pseudochirops_corinnae_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Pseudochirops_cupreus-GCA_028627135.1_Pseudochirops_cupreus_HiC_genomic.fna $i $chr_id &
./genome_chain_toga.sh Thylacinus_cynocephalus-GCA_007646695.3_UniMelb_ThyCyn2.0_hybrid_assembly_genomic.fna $i $chr_id &
./genome_chain_toga.sh Trichosurus_vulpecula-GCA_011100635.1_mTriVul1.pri_genomic.fna $i $chr_id &
./genome_chain_toga.sh Vombatus_ursinus-GCA_028626985.1_vu-2k_genomic.fna $i $chr_id &

# Wait for all parallel jobs to finish
wait

# Continue with the next steps
echo "All jobs for $i completed."
