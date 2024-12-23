wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/635/505/GCF_902635505.1_mSarHar1.11/GCF_902635505.1_mSarHar1.11_genomic.gtf.gz
gunzip GCF_902635505.1_mSarHar1.11_genomic.gtf.gz
awk '$3=="CDS"' GCF_902635505.1_mSarHar1.11_genomic.gtf > GCF_902635505.1_mSarHar1.11_genomic.gtf.cds.gtf
