#!/bin/bash
# This script
# -- dowloads to raw_data dir:
# 1) reference E. coli genome in fasta format (genomic.fna.gz);
# 2) reference E. coli genome annotation in GFF3 format (genomic.gff.gz);
# 3) files with MD5 checksums (for archives and uncompressed files as well);
# 4) raw shotgun Illumina PE reads (_R1, _R2) of ampicillin-resistant E. coli stain in fastq format (fastq.gz).
# -- checks MD5 sums for:
# 1) reference and annotation (md5checksums.txt);
# 2) raw reads _R1 and _R2 (md5checksums_amp_res_reads.txt).
# -- writes to logs dir:
# 1) data_downloading.log
# 2) md5_check.log


# exit on error, undefined variables, and pipe failures
set -e
set -u
set -o pipefail

##-- DOWNLOAD --

echo "Downloading reference genome, annotation file and md5checksums..."
wget -c -P raw_data -a logs/data_downloading.log \
	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz \
	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz \
	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/md5checksums.txt \
	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/uncompressed_checksums.txt

# downlad raw reads from figshare
# pretend to be a browser Mozilla/5.0 to pass figshare block
# wget writes downloading progress messages to stderr 2>

echo "Downloading reads _R1.fastq.gz ..."
wget -c --user-agent="Mozilla/5.0" \
  -O raw_data/amp_res_R1.fastq.gz \
  https://figshare.com/ndownloader/files/23769689 \
  -a logs/data_downloading.log

echo "Downloading reads _R2.fastq.gz ..."
wget -c --user-agent="Mozilla/5.0" \
  -O raw_data/amp_res_R2.fastq.gz \
  https://figshare.com/ndownloader/files/23769692 \
  -a logs/data_downloading.log

##--CHECK MD5 CHECKSUMS --
# stdout and stderr directed to log file
# grep otherwise the script will exit cause FAILED checks on un-downloaded files related to reference
echo "Checking MD5 checksums for reference genome and genome annotation"
cd raw_data
grep -E "ASM584v2_genomic\.(fna|gff)\.gz" md5checksums.txt | md5sum -c > ../logs/md5_check.log 2>&1 
cd ..

# create md5checksums file with hashes from figshare info and file names
# field separator for md5checksums files is 2 spaces

cat > raw_data/md5checksums_amp_res_reads.txt <<EOF 
181588f2e9196486479381bf3a0611e9  amp_res_R1.fastq.gz
2a768a51991ac35baf8847c9be2a51f3  amp_res_R2.fastq.gz
EOF

echo "" >> logs/md5_check.log
echo "Created MD5 checksums file for ampicillin-resistant reads:" >> logs/md5_check.log
cat raw_data/md5checksums_amp_res_reads.txt >> logs/md5_check.log

echo "Checking MD5 checksums for _R1 and _R2.fastq.gz read files"
echo "" >> logs/md5_check.log
cd raw_data
md5sum -c md5checksums_amp_res_reads.txt >> ../logs/md5_check.log 2>&1
cd ..

echo "Download and MD5 verification completed successfully!"
