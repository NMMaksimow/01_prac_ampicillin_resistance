# LabJournal
# Nikita Maksimov
# Project: Identification of ampicillin-resistance mechanism in *E. coli*
# Goal: Identify mutations conferring ampicillin resistance by comparing resistant strain to reference genome
----

#### :calendar: 10 November 2025
### Project seting up project directory, downloading data, checking MD5 checksums and inspecing files
<details>
<summary>Basic commands to activate mamba environment, set up directory tree for the project, inititate git and create .gitignore</summary>

Activate environment:
```bash
mamba activate bioinf
```
Create project directory tree:
```bash
mkdir -p 01_prac_ampicillin_resistance/{raw_data,logs,scripts,results}
cd 01_prac_ampicillin_resistance/
```
Initiate Git:
```bash
git init
```
Create .gitignore using Here Document and End Of File (EOF text text EOF):
```bash
cat > .gitignore <<EOF
# Raw data - too large for git
raw_data/*.gz
raw_data/*.fna
raw_data/*.fasta
raw_data/*.fastq
raw_data/*.bam
raw_data/*.sam

# Results - generated during analysis
results/

# Python
__pycache__/
*.pyc

# R
.Rhistory
.RData

# Editor
*.swp
*~
.DS_Store
EOF
```
</details>

After set up I made an initial Git commit.

Programmatic download of reference genome, its annotation and MD5 checksums were straightforward with `wget` command. Download of read files using wget were a bit dodgy because of access to figshare (403 Forbidden).

First I thought it's a connection issue and used these commands for interactive check it.

```bash
wget --spider https://figshare.com/ndownloader/files/23769689
wget --spider https://figshare.com/ndownloader/files/23769692
```

Then I figured it out that Figshare blocks programmatic access and work around it using `--user-agent` option for wget (in script below). To check MD5 checksum I used md5checksums.txt (for reference genome and annotation) and created md5checksums_amp_res_reads.txt from MD5 on Figshare (for reads). I used `grep` (with -E flag for Extended RegEx) to filter only downloaded files from others present on ncbi (use backslash to escape the dot: \. matches a literal period, not "any character" as in regex). 

Full script: [`scripts/01_downloading_raw_data_md5_check.sh`](scripts/01_downloading_raw_data_md5_check.sh)

To run the script I set up execution rights and use standard command:
```bash
chmod +x scripts/01_downloading_raw_data_md5_check.sh
bash scripts/01_downloading_raw_data_md5_check.sh
```
<details>
<summary>The console output</summary>
(bioinf) ➜  01_prac_ampicillin_resistance git:(master) bash scripts/01_downloading_raw_data_md5_check.sh
Downloading reference genome, annotation file and md5checksums...
Downloading reads _R1.fastq.gz ...
Downloading reads _R2.fastq.gz ...
Checking MD5 checksums for reference genome and genome annotation
Checking MD5 checksums for _R1 and _R2.fastq.gz read files
Download and MD5 verification completed successfully!
</details>

I performed Git commit after downloading data and its integrity.

I adore tree command to check the project directory structure, it's very visual. I'll try not to use it too much.
```bash
tree .
```

```
.
├── 01_lab_journal_nm.txt
├── README.md
├── logs
│   ├── data_downloading.log
│   └── md5_check.log
├── raw_data
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz
│   ├── GCF_000005845.2_ASM584v2_genomic.gff.gz
│   ├── amp_res_R1.fastq.gz
│   ├── amp_res_R2.fastq.gz
│   ├── md5checksums.txt
│   ├── md5checksums_amp_res_reads.txt
│   └── uncompressed_checksums.txt
├── results
└── scripts
    └── 01_downloading_raw_data_md5_check.sh

5 directories, 12 files
```

I used zcat and head commands to verify reference genome fna.gz data format. The FASTA has header line starting wiht > symbol and DNA sequence as expected :white_check_mark:

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz | head -20
```

<details>
<summary>Console output: FASTA header and 19 more lines</summary>
>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG
GTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGAC
AGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT
AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGG
TAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCG
ATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTG
GCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT
GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAA
AACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAA
ATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCT
GGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTG
ATCACATGGTGCTGATGGCAGGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGAC
TACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTG
CGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCG
CTAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATACCGGAAATCCT
CAAGCACCAGGTACGCTCATTGGTGCCAGCCGTGATGAAGACGAATTACCGGTCAAGGGCATTTCCAATCTGAATAACAT
GGCAATGTTCAGCGTTTCTGGTCCGGGGATGAAAGGGATGGTCGGCATGGCGGCGCGCGTCTTTGCAGCGATGTCACGCG
CCCGTATTTCCGTGGTGCTGATTACGCAATCATCTTCCGAATACAGCATCAGTTTCTGCGTTCCACAAAGCGACTGTGTG
CGAGCTGAACGGGCAATGCAGGAAGAGTTCTACCTGGAACTGAAAGAAGGCTTACTGGAGCCGCTGGCAGTGACGGAACG
</details>

I used zcat, grep, wc commands to check how many contigs are in the reference genome. Wordcount command was used with option lines (-l). There's only one contig -- bacterial circular chromosome as expected.:white_check_mark:

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz | grep -E "^>" | wc -l

# Output: 1
```

To calculate reference genome size in basepairs I used pipe of zcat, grep (to filter the header of), transliterate with delete option (to get rid of newline character) and wordcount with character option. The reference genome size is 4,641,652 bp (4,6 Mb).
The sanity check is passed. *E. coli* genomes size vary: 4,5 — 5.5 Mb, depending on strain. K-12 is a lab strain that might expirienced size reduction being passed in lab environmnet for miriads of generations. :white_check_mark:

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz | grep -v "^>" | tr -d '\n' | wc -c

# Output: 4641652
```

Again, I used zcat and head pipe to verify genome annotation in GFF3 format. All looks as expected. :white_check_mark:

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.gff.gz | head -20
```

<details>
<summary>Console output: genofe annotation GFF-3 head</summary>

##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
#!genome-build ASM584v2
#!genome-build-accession NCBI_Assembly:GCF_000005845.2
##sequence-region NC_000913.3 1 4641652
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=511145
NC_000913.3     RefSeq  region  1       4641652 .       +       .       ID=NC_000913.3:1..4641652;Dbxref=taxon:511145;Is_circular=true;Name=ANONYMOUS;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=K-12;substrain=MG1655
NC_000913.3     RefSeq  gene    190     255     .       +       .       ID=gene-b0001;Dbxref=ASAP:ABE-0000006,ECOCYC:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001;locus_tag=b0001
NC_000913.3     RefSeq  CDS     190     255     .       +       0       ID=cds-NP_414542.1;Parent=gene-b0001;Dbxref=UniProtKB/Swiss-Prot:P0AD86,GenBank:NP_414542.1,ASAP:ABE-0000006,ECOCYC:EG11277,GeneID:944742;Name=NP_414542.1;gbkey=CDS;gene=thrL;locus_tag=b0001;product=thr operon leader peptide;protein_id=NP_414542.1;transl_table=11
NC_000913.3     RefSeq  gene    337     2799    .       +       .       ID=gene-b0002;Dbxref=ASAP:ABE-0000008,ECOCYC:EG10998,GeneID:945803;Name=thrA;gbkey=Gene;gene=thrA;gene_biotype=protein_coding;gene_synonym=ECK0002,Hs,thrA1,thrA2,thrD;locus_tag=b0002
</details>

GFF3 starts with header lines (marked with #) and then the body consists tab-separated fields:

`seqid | source | type | start | end | score | strand | phase | attributes`

- `seqid` — Chromosome/contig ID
- `source` — Annotation source (e.g., RefSeq)
- `type` — Feature type (gene, CDS, exon, etc.)
- `start` — Start position (1-based)
- `end` — End position
- `score` — Score (`.` if not applicable)
- `strand` — `+` or `-`
- `phase` — Reading frame for CDS (0, 1, 2)
- `attributes` — Semicolon-separated key=value pairs


I used zcat, grep and wordcount to count all entries in GFF genome annotation. I used invert match option (-v "^#") to remove header lines in annotation before counting.

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.gff.gz | grep -v "^#" | wc -l

# Output: 9523
```

<details>
<summary>Without grep filtering I got total number of lines in GFF annotation including header lines.</summary>

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.gff.gz | wc -l

# Output: 9531
```
</details>

To further explore GFF annotation I used zcat, grep, cut, sort, wc pipe to examine the unique feature types (the third column). I used option field (-f3) to specify third column for cut command. There are 11 unique features in reference annotation.

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.gff.gz| grep -v "^#" | cut -f3 | sort -u | wc -l

# Output: 11
```

To list features in annotation I combined zcat, grep, cut, sort and uniq with count flag (-c).

```bash
zcat raw_data/GCF_000005845.2_ASM584v2_genomic.gff.gz | grep -v "^#" | cut -f3 | sort | uniq -c

# Output:
#   4340 CDS
#    216 exon
#   4506 gene
#     50 mobile_genetic_element
#    108 ncRNA
#      1 origin_of_replication
#    145 pseudogene
#     22 rRNA
#      1 region
#     48 sequence_feature
#     86 tRNA
```

To count R1 and R2 reads in FASTQ archives I used pipe of zcat, wc, awk commands. I divided the wordcount output by 4. Each read is represented by 4 lines: header, sequence, plus placeholder, quality control. There are 455,876 of forward (R1) and reverse (R2) reads. All reads are paired.

```bash
zcat raw_data/amp_res_R1.fastq.gz | wc -l | awk '{print $1/4}'

# Output: 455876
zcat raw_data/amp_res_R2.fastq.gz | wc -l | awk '{print $1/4}'

# Output: 455876
```

I used seqkit to collect summary statistics on FASTA and FASTQ files. I installed seqkit to bioinf environment using mamba, specifying bioconda channel with option -c. 

```bash
mamba install -c bioconda seqkit
```

I run seqkit on FASTA genome reference. The seqkit output corroborated manual check. The FASTA file has only one ≈4,6 Mb contig (circular bacterial chromosome).

```bash
seqkit stats raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz

# Output:
# file                                              format  type  num_seqs    sum_len    min_len    avg_len    max_len
# raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz  FASTA   DNA          1  4,641,652  4,641,652  4,641,652  4,641,652
```

Running seqkit on read FASTQ files also verified results of read counting with in-built commands. All reads have equal length of 101 bp.

```bash
seqkit stats raw_data/amp_res_R*.fastq.gz

# file                          format  type  num_seqs     sum_len  min_len  avg_len  max_len
# raw_data/amp_res_R1.fastq.gz  FASTQ   DNA    455,876  46,043,476      101      101      101
# raw_data/amp_res_R2.fastq.gz  FASTQ   DNA    455,876  46,043,476      101      101      101
```

### Quality control using FastQC
<details>
<summary> I installed FastQC using mamba package manager.</summary>

```bash
mamba install -c bioconda fastqc
fastqc --version

# Output: FastQC v0.12.1
```
</details>

I performed quality control using FastQC with the following command. I redirected both stdout and stderr to log file.
```bash
mkdir -p results/fastqc
fastqc raw_data/*.fastq.gz -o results/fastqc &> logs/fastqc.log
```
There is a typical to Illumina reads quality drop at the end of the reads.

### Quality-based trimming and post-trimming quality control with FastQC
#### fastp trimming
<details>
<summary>I used fastp for trimming at the beginning, because I had expirience working with it before.</summary>

I installed fastp using mamba and created directory for trimmed files.
```bash
mamba install -c bioconda fastqc

fastp --version

# Output: fastp 1.0.1

mkdir -p results/fastp_trimmed
```
</details>

After inspecting per base sequence quality plots I run fastp using the following flags and options. It was pretty relaxed trimming with >=20 quality threshold (option `-q`, `--qualified_quality_phred`), 25% of bases are allowed to be unqualified (option `-u`, `--unqualified_percent_limit`). All reads with >5 N bases were discarded (option `-n`, `--n_base_limit`). All reads shorter than 50 bp were discarded (option `-l`, `--length_required`). I used flag `--detect_adapter_for_pe`. I generated JSON and HTML repots. I redirected stderr to the same log file where stdout was written.

```bash
fastp \
  -i raw_data/amp_res_R1.fastq.gz \
  -I raw_data/amp_res_R2.fastq.gz \
  -o results/fastp_trimmed/amp_res_R1_trimmed.fastq.gz \
  -O results/fastp_trimmed/amp_res_R2_trimmed.fastq.gz \
  -q 20 -u 25 -n 5 -l 50 \
  --detect_adapter_for_pe \
  --html results/fastp_trimmed/fastp_amp_res.html \
  --json results/fastp_trimmed/fastp_amp_res.json \
  > logs/fastp_trimming.log 2>&1
```

After trimming with fastp I run fastqc on trimmed files using the following command.
```bash
mkdir -p results/fastqc_fastp_trimmed
fastqc results/fastp_trimmed/*.fastq.gz -o results/fastqc_fastp_trimmed &> logs/fastqc_fastp_trimmed.log
```

The number of trimmed sequences is equal for _R1 and _R2 (all paired): 380,531 reads passed the filter.

#### trimmomatic relaxed trimming

<details>
<summary>I installed trimmomatic in the bioinf environment and created directory for output files.</summary>

```bash
mamba install -c bioconda Trimmomatic
trimmomatic -version

# Output: 0.40

mkdir -p results/trimmomatic_trimmed
```
</details>

To run trimmomatic I used flags PE (for paired end reads) and `-phred33` to specify QS encoding, because FastQC reports specify Illumina machine: "Encoding	Sanger / Illumina 1.9". Trimmomatic output 4 files: unpaired reads (_1U and _2U) and paired reads that passed filter (_1P and _2P). LEADING:20 — trim bases from the start of a read if quality < 20. TRAILING:20 — trim bases from the end of a read if quality < 20. SLIDINGWINDOW:10:20 — scan with a 10-bp window, cut when average quality drops below 20. MINLEN:20 — discard reads shorter than 20 bp after trimming.

```bash
trimmomatic PE -threads 4 -phred33 \
  raw_data/amp_res_R1.fastq.gz raw_data/amp_res_R2.fastq.gz \
  -baseout results/trimmomatic_trimmed/amp_res.fastq.gz \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20 \
  > logs/trimmomatic_trimming.log 2>&1
 ```

Adter trimming using trimmomatic I run FastQC on trimmed files.
```bash
mkdir -p results/fastqc_trimmomatic_trimmed
fastqc results/trimmomatic_trimmed/*P.fastq.gz -o results/fastqc_trimmomatic_trimmed &> logs/fastqc_trimmomatic_trimmed.log
```
 
#### trimmomatic strict trimming
For strict trimming I rise quality thresholds to 30. 
```bash
mkdir -p results/trimmomatic_strictly_trimmed

trimmomatic PE -threads 4 -phred33 \
  raw_data/amp_res_R1.fastq.gz raw_data/amp_res_R2.fastq.gz \
  -baseout results/trimmomatic_strictly_trimmed/amp_res_strictly.fastq.gz \
  LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:20 \
  > logs/trimmomatic_strict_trimming.log 2>&1
```

FastQC on paired fastq files after strict trimming with trimmomatic:
```bash
mkdir -p results/fastqc_trimmomatic_strictly_trimmed
fastqc results/trimmomatic_strictly_trimmed/*P.fastq.gz -o results/fastqc_trimmomatic_strictly_trimmed &> logs/fastqc_trimmomatic_strictrly_trimmed.log
```

### Alignment and mapping
bwa already exists in the environment, so there's no need for installation

Creating index for reference genome using bwa index command:
```bash
bwa index raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz &> logs/bwa_index.log
```

NB! Resulting BW index files (suffixes .amb, .ann, .bwt, .pac, .sa) are located in raw_data directory, so that other programs can find them.

├── raw_data
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.amb
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.ann
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.bwt
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.pac
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.sa
│   ├── GCF_000005845.2_ASM584v2_genomic.gff.gz
│   ├── amp_res_R1.fastq.gz
│   ├── amp_res_R2.fastq.gz
│   ├── md5checksums.txt
│   ├── md5checksums_amp_res_reads.txt
│   └── uncompressed_checksums.txt

Aligning and mapping trimmed R1 and R2 reads to reference genome:

```bash
mkdir results/bwa_alignments
```

Three alignments SAM files will be generated based on three different trimming strategies.
# trimmomatic trimming
```bash
bwa mem raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz \
	results/trimmomatic_trimmed/amp_res_1P.fastq.gz \
	results/trimmomatic_trimmed/amp_res_2P.fastq.gz \
	> results/bwa_alignments/aln_trimmomatic_trimmed.sam 2> logs/bwa_mem_trimmomatic.log
```
# strict trimmomatic trimming
```bash	
bwa mem raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz \
	results/trimmomatic_strictly_trimmed/amp_res_strictly_1P.fastq.gz \
	results/trimmomatic_strictly_trimmed/amp_res_strictly_2P.fastq.gz \
	> results/bwa_alignments/aln_trimmomatic_strictly_trimmed.sam 2> logs/bwa_mem_trimmomatic_strictly.log \
```
# fastp trimming
```bash
bwa mem raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz \
	results/fastp_trimmed/amp_res_R1_trimmed.fastq.gz \
	results/fastp_trimmed/amp_res_R2_trimmed.fastq.gz \
	> results/bwa_alignments/aln_fastp.sam 2> logs/bwa_mem_fastp.log 
```

├── results
│   ├── bwa_alignments
│   │   ├── aln_fastp.sam
│   │   ├── aln_trimmomatic_strictly_trimmed.sam
│   │   └── aln_trimmomatic_trimmed.sam

Inspect SAM files:
```bash
less -S results/bwa_alignments/aln_trimmomatic_strictly_trimmed.sam
```
For loop for: converting SAM to binary alignment map file (BAM), sorting BAM (add suffix _sorted.bam) and building index on sorted files (.bai)
```bash
for sam in results/bwa_alignments/*.sam; do
    # basename 
    base=$(basename "$sam" .sam)

    echo "Processing $base ..."

    # convert SAM format alignment to binary version (BAM format)
    samtools view -b "$sam" > results/bwa_alignments/${base}.bam 2>> logs/samtools_samtobam_sorting_indexing.log

    # sorting BAM files before indexing
    samtools sort results/bwa_alignments/${base}.bam \
        -o results/bwa_alignments/${base}_sorted.bam 2>> logs/samtools_samtobam_sorting_indexing.log

    # indexing BAM
    samtools index results/bwa_alignments/${base}_sorted.bam 2>> logs/samtools_samtobam_sorting_indexing.log

    echo "$base done." >> logs/samtools_samtobam_sorting_indexing.log
done
```

IGV need uncompressed reference genome in FASTA format so I extracted reference from the archive using gunzip.
Flag -k preserves original fna.gz archive.
```bash
gunzip -k raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

I decided to check md5 checksum. But the command I used earlier failed.
```bash
md5sum -c raw_data/uncompressed_checksums.txt
```

Checksums for uncompressed files from ncbi is stored in 4 column format.
(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ cat raw_data/uncompressed_checksums.txt
#file   md5sum  crc32   size
./GCF_000005845.2_ASM584v2_ani_contam_ranges.tsv        0feb9ab0ef81aed6d768628d50b119fd        0c012c82        46803
./GCF_000005845.2_ASM584v2_ani_report.txt       e3cc1640429bfbb576a1328b0b371fff        bb0bd3b6        5812
./GCF_000005845.2_ASM584v2_assembly_report.txt  5d767338dae878d0d547ca5edb138df3        28f934b1        1207
./GCF_000005845.2_ASM584v2_cds_from_genomic.fna 113ad98b89cb697e22aae117a1ae1a89        b0a79c52        5043627
./GCF_000005845.2_ASM584v2_fcs_report.txt       9bd2c1f0b8d87e8b1b7215551c099c17        924e232e        667
./GCF_000005845.2_ASM584v2_genomic.fna  92c997bcd88e983ffdb21b2712ed3736        f36bbf76        4699745
./GCF_000005845.2_ASM584v2_genomic.gbff 80edc82b4a5d19198f1eee5a3848fdac        c80a55e4        11882100
./GCF_000005845.2_ASM584v2_genomic.gff  a7cef966c981b35c5bbf30c3df923312        f5f907ce        2467646
./GCF_000005845.2_ASM584v2_genomic.gtf  27882704f8ef02338681dcdd2ee6d26b        4bbbd7df        6774966
./GCF_000005845.2_ASM584v2_protein.faa  4c19a285fd7a4e75ce34a466fcf0607f        25917a05        1750878
./GCF_000005845.2_ASM584v2_rna_from_genomic.fna fd55823557f1ecb83c3d8939a0e88aea        0796f039        97814
./GCF_000005845.2_ASM584v2_translated_cds.faa   1d55819ca9d1172f528ff26ccb9687a0        7b546c56        2325398

Commmand to quick check of md5checksums on uncompressed files.
awk changes the order of hashes and file names to feed md5sum
tail -n +2 command start printin from the second line and reads input up to the end of file (need to omit headers from 4-column ncbi checksum file)
```bash
cd raw_data
awk '{print $2, $1}' uncompressed_checksums.txt | tail -n +2 | md5sum -c
cd ..
```

(bioinf) ➜  raw_data git:(master) ✗ awk '{print $2, $1}' uncompressed_checksums.txt | tail -n +2 | md5sum -c
md5sum: ./GCF_000005845.2_ASM584v2_ani_contam_ranges.tsv: No such file or directory
./GCF_000005845.2_ASM584v2_ani_contam_ranges.tsv: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_ani_report.txt: No such file or directory
./GCF_000005845.2_ASM584v2_ani_report.txt: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_assembly_report.txt: No such file or directory
./GCF_000005845.2_ASM584v2_assembly_report.txt: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_cds_from_genomic.fna: No such file or directory
./GCF_000005845.2_ASM584v2_cds_from_genomic.fna: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_fcs_report.txt: No such file or directory
./GCF_000005845.2_ASM584v2_fcs_report.txt: FAILED open or read
./GCF_000005845.2_ASM584v2_genomic.fna: OK
md5sum: ./GCF_000005845.2_ASM584v2_genomic.gbff: No such file or directory
./GCF_000005845.2_ASM584v2_genomic.gbff: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_genomic.gff: No such file or directory
./GCF_000005845.2_ASM584v2_genomic.gff: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_genomic.gtf: No such file or directory
./GCF_000005845.2_ASM584v2_genomic.gtf: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_protein.faa: No such file or directory
./GCF_000005845.2_ASM584v2_protein.faa: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_rna_from_genomic.fna: No such file or directory
./GCF_000005845.2_ASM584v2_rna_from_genomic.fna: FAILED open or read
md5sum: ./GCF_000005845.2_ASM584v2_translated_cds.faa: No such file or directory
./GCF_000005845.2_ASM584v2_translated_cds.faa: FAILED open or read
md5sum: WARNING: 11 listed files could not be read

Inspect in IGV:
1) genome reference in .fna,
2) archived annotation gff.gz and 
3) alignment map _sorted.bam with index _sorted.bam.bai in the same directory

Some summary on alignment map files.
```bash
mkdir -p results/alignment_stats
```

```bash
for bam in results/bwa_alignments/*_sorted.bam; do
    base=$(basename "$bam" _sorted.bam)
    samtools flagstat "$bam" > results/alignment_stats/${base}_flagstat.txt
    samtools idxstats "$bam" > results/alignment_stats/${base}_idxstats.txt
done
```

##############################################
#### Script to collect alignment stats #######
##############################################
#!/usr/bin/env bash
set -euo pipefail

# === PATHS ===
ALIGN_DIR="results/bwa_alignments"
STATS_DIR="results/alignment_stats"
LOG_DIR="logs"
OUTFILE="${STATS_DIR}/alignment_summary.tsv"

mkdir -p "$STATS_DIR" "$LOG_DIR"

echo -e "sample\ttotal_reads\tmapped_reads\tproperly_paired\tpercent_mapped\tpercent_proper\tmean_depth\tavg_MAPQ" > "$OUTFILE"

# === MAIN LOOP ===
for bam in ${ALIGN_DIR}/*_sorted.bam; do
    base=$(basename "$bam" _sorted.bam)
    flagfile="${STATS_DIR}/${base}_flagstat.txt"
    idxfile="${STATS_DIR}/${base}_idxstats.txt"

    echo "Processing ${base}..."

    # total reads
    total=$(grep "in total" "$flagfile" | awk '{print $1}')
    # mapped reads
    mapped=$(grep " mapped (" "$flagfile" | head -n1 | awk '{print $1}')
    # properly paired reads
    proper=$(grep "properly paired" "$flagfile" | awk '{print $1}')
    # % mapped
    pmapped=$(awk -v m=$mapped -v t=$total 'BEGIN{printf "%.2f", (m/t)*100}')
    # % proper
    pproper=$(awk -v p=$proper -v t=$total 'BEGIN{printf "%.2f", (p/t)*100}')

    # mean depth (все позиции, включая нулевые)
    depth=$(samtools depth -a "$bam" | awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}')

    # среднее значение MAPQ
    mapq=$(samtools view "$bam" | awk '{m[$5]++} END{for (k in m){sum+=k*m[k]; n+=m[k]} if(n>0) printf "%.2f", sum/n; else print 0}')

    echo -e "${base}\t${total}\t${mapped}\t${proper}\t${pmapped}\t${pproper}\t${depth}\t${mapq}" >> "$OUTFILE"
done

echo "Summary written to: $OUTFILE"

##############################################
##############################################


## ---- 06 Variant calling ----
```bash
mamba install -c bioconda -c conda-forge varscan
```
```bash
varscan --help
```

```bash
mkdir -p results/variant_calling
```
Create mpileup files to count how many reads have a mutation at the same position.
mpileup files are required for SNP calling by varscan

```bash
for bam in results/bwa_alignments/*_sorted.bam; do
    base=$(basename "$bam" _sorted.bam)
    echo "Calling variants for $base..."

    samtools mpileup \
        -f raw_data/GCF_000005845.2_ASM584v2_genomic.fna \
        -o results/variant_calling/${base}.mpileup \
        "$bam" \
        2>> logs/samtools_mpileup.log
done
```
Extracts on some of samtools mpileup options from samtools mpileup --help:
  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)
  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
 
SNP calling
```bash
for mpileup in results/variant_calling/*.mpileup; do
    base=$(basename "$mpileup" .mpileup)
    echo "Identifying SNPs from $base ..."
    
    varscan mpileup2snp "$mpileup" \
        --min-var-freq 0.7 \
        --variants \
        --output-vcf 1 \
        > results/variant_calling/varscan_res_${base}.vcf \
        2>> logs/varscan_mpileup2snp.log
done
```

VCF (variant call format) is a tabular format and it could be digested by IGV for visualization.
(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ less -S results/variant_calling/varscan_res_aln_trimmomatic_trimmed.vcf | grep -v "^##" | cat
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1
NC_000913.3     93043   .       C       G       .       PASS    ADP=17;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:93:18:17:0:17:100%:4.2852E-10:0:36:0:0:7:10
NC_000913.3     482698  .       T       A       .       PASS    ADP=16;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:87:16:16:0:16:100%:1.6637E-9:0:45:0:0:7:9
NC_000913.3     852762  .       A       G       .       PASS    ADP=14;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:76:14:14:0:14:100%:2.4927E-8:0:36:0:0:8:6
NC_000913.3     1905761 .       G       A       .       PASS    ADP=13;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:70:13:13:0:13:100%:9.6148E-8:0:44:0:0:11:2
NC_000913.3     3535147 .       A       C       .       PASS    ADP=17;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:93:17:17:0:17:100%:4.2852E-10:0:36:0:0:10:7
NC_000913.3     4390754 .       G       T       .       PASS    ADP=15;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:81:16:15:0:15:100%:6.4467E-9:0:36:0:0:8:7
(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ less -S results/variant_calling/varscan_res_aln_trimmomatic_strictly_trimmed.vcf | grep -v "^##" | cat
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1
NC_000913.3     93043   .       C       G       .       PASS    ADP=9;WT=0;HET=0;HOM=1;NC=0     GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:46:9:9:0:9:100%:2.0568E-5:0:36:0:0:3:6
NC_000913.3     482698  .       T       A       .       PASS    ADP=13;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:70:13:13:0:13:100%:9.6148E-8:0:43:0:0:5:8
NC_000913.3     852762  .       A       G       .       PASS    ADP=9;WT=0;HET=0;HOM=1;NC=0     GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:46:9:9:0:9:100%:2.0568E-5:0:38:0:0:6:3
NC_000913.3     1905761 .       G       A       .       PASS    ADP=11;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:58:11:11:0:11:100%:1.4176E-6:0:43:0:0:9:2
NC_000913.3     3535147 .       A       C       .       PASS    ADP=9;WT=0;HET=0;HOM=1;NC=0     GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:46:9:9:0:9:100%:2.0568E-5:0:32:0:0:4:5
NC_000913.3     4390754 .       G       T       .       PASS    ADP=12;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:64:12:12:0:12:100%:3.698E-7:0:37:0:0:5:7
(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ less -S results/variant_calling/varscan_res_aln_fastp.vcf | grep -v "^##" | cat
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1
NC_000913.3     93043   .       C       G       .       PASS    ADP=13;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:70:13:13:0:13:100%:9.6148E-8:0:36:0:0:4:9
NC_000913.3     482698  .       T       A       .       PASS    ADP=15;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:81:15:15:0:15:100%:6.4467E-9:0:45:0:0:6:9
NC_000913.3     852762  .       A       G       .       PASS    ADP=13;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:70:13:13:0:13:100%:9.6148E-8:0:36:0:0:8:5
NC_000913.3     1905761 .       G       A       .       PASS    ADP=12;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:64:12:12:0:12:100%:3.698E-7:0:45:0:0:10:2
NC_000913.3     3535147 .       A       C       .       PASS    ADP=10;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:52:10:10:0:10:100%:5.4125E-6:0:39:0:0:5:5
NC_000913.3     4390754 .       G       T       .       PASS    ADP=16;WT=0;HET=0;HOM=1;NC=0    GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    1/1:87:16:16:0:16:100%:1.6637E-9:0:34:0:0:8:8


## ---- 07 Automatic SNP annotation and variant effect prediction ----

```bash
# downlad snpEff for automatica annotation
mamba install -c bioconda snpEff

# create a directory with K12 reference strain database
mkdir -p results/snpEff/data/k12

# create snpEff.config file in variant_calling directory
cat > results/snpEff/snpEff.config <<EOF
k12.genome : ecoli_K12
EOF

# download GBFF file with annotation and reference sequence
wget -c -P results/snpEff/data/k12 -a logs/gbff_download.log \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz

# unzip file with renaming to database folder
gzip -d results/snpEff/data/k12/*.gz
mv results/snpEff/data/k12/*.gbff results/snpEff/data/k12/genes.gbk

# create database
snpEff build -genbank -v k12 -c results/snpEff/snpEff.config &> logs/snpEff_build.log

# for loop for annotation of three different VCF files
for vcf in results/variant_calling/*.vcf; do
    base=$(basename "$vcf" .vcf)
    echo "Annotating $base.vcf ..."
    
    snpEff ann -c results/snpEff/snpEff.config k12 "$vcf" \
        > results/variant_calling/${base}_annotated.vcf \
        2>> logs/snpEff_annotation.log
done
```

```bash
grep "ANN=" results/variant_calling/varscan_res_aln_fastp_annotated.vcf \
  | sed 's/.*ANN=//' | cut -d'|' -f4 | sort | uniq -c
```

(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ grep "ANN=" results/variant_calling/varscan_res_aln_fastp_annotated.vcf \
  | sed 's/.*ANN=//' | cut -d'|' -f4 | sort | uniq -c

      1 acrB
      1 envZ
      1 ftsI
      1 glnH
      1 mntP
      1 rsgA
	  
################################################################
################ 11/11/2025 ####################################
################################################################

# ============================================================
# ============================================================

#!/usr/bin/env bash
# ============================================================
# 00_check_software_versions.sh
# Logs versions of all bioinformatics tools used in the pipeline
# ============================================================

set -euo pipefail
trap 'echo "Error at line $LINENO in $0" >&2' ERR

# --- Setup ---
LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/software_versions.log"

echo "Environment freeze for E. coli ampicillin resistance project - $(date)" > "$LOG_FILE"

# --- Helper function ---
log_version() {
    local name="$1"
    local cmd="$2"
    echo "${name} version:" >> "$LOG_FILE"
    if eval "$cmd" >> "$LOG_FILE" 2>&1; then
        echo "" >> "$LOG_FILE"
    else
        echo "  [Warning] ${name} not found or failed to report version." >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"
    fi
}

# --- Core Tools ---
log_version "Python" "python --version"
log_version "mamba" "mamba --version"
log_version "bwa" "bwa 2>&1 | grep -m1 'Program'"
log_version "samtools" "samtools --version | head -n1"
log_version "fastqc" "fastqc --version"
log_version "fastp" "fastp --version"
log_version "Trimmomatic" "trimmomatic -version"
log_version "seqkit" "seqkit version 2>&1 | head -n1"
log_version "varscan" "varscan --help | grep -m1 'VarScan'"
log_version "snpEff" "snpEff -version 2>&1 | head -n1"

# --- Optional dependencies ---
log_version "git" "git --version"
log_version "wget" "wget --version | head -n1"
log_version "gzip" "gzip --version | head -n1"

echo "Software version logging complete."
echo "Log saved to: ${LOG_FILE}"
# ============================================================

#######################################

(base) ➜  01_prac_ampicillin_resistance git:(master) ✗ tree .
.
├── 01_lab_journal_nm.txt
├── README.md
├── logs
│   ├── bwa_index.log
│   ├── bwa_mem_fastp.log
│   ├── bwa_mem_trimmomatic.log
│   ├── bwa_mem_trimmomatic_strictly.log
│   ├── data_downloading.log
│   ├── fastp_trimming.log
│   ├── fastqc.log
│   ├── fastqc_fastp_trimmed.log
│   ├── fastqc_trimmomatic_strictrly_trimmed.log
│   ├── fastqc_trimmomatic_trimmed.log
│   ├── gbff_download.log
│   ├── md5_check.log
│   ├── samtools_mpileup.log
│   ├── samtools_samtobam_sorting_indexing.log
│   ├── snpEff_annotation.log
│   ├── snpEff_build.log
│   ├── software_versions.log
│   ├── trimmomatic_strict_trimming.log
│   ├── trimmomatic_trimming.log
│   └── varscan_mpileup2snp.log
├── raw_data
│   ├── GCF_000005845.2_ASM584v2_genomic.fna
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.fai
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.amb
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.ann
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.bwt
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.pac
│   ├── GCF_000005845.2_ASM584v2_genomic.fna.gz.sa
│   ├── GCF_000005845.2_ASM584v2_genomic.gff.gz
│   ├── amp_res_R1.fastq.gz
│   ├── amp_res_R2.fastq.gz
│   ├── md5checksums.txt
│   ├── md5checksums_amp_res_reads.txt
│   └── uncompressed_checksums.txt
├── results
│   ├── alignment_stats
│   │   ├── alignment_summary.tsv
│   │   ├── aln_fastp_flagstat.txt
│   │   ├── aln_fastp_idxstats.txt
│   │   ├── aln_trimmomatic_strictly_trimmed_flagstat.txt
│   │   ├── aln_trimmomatic_strictly_trimmed_idxstats.txt
│   │   ├── aln_trimmomatic_trimmed_flagstat.txt
│   │   └── aln_trimmomatic_trimmed_idxstats.txt
│   ├── bwa_alignments
│   │   ├── aln_fastp.bam
│   │   ├── aln_fastp.sam
│   │   ├── aln_fastp_sorted.bam
│   │   ├── aln_fastp_sorted.bam.bai
│   │   ├── aln_trimmomatic_strictly_trimmed.bam
│   │   ├── aln_trimmomatic_strictly_trimmed.sam
│   │   ├── aln_trimmomatic_strictly_trimmed_sorted.bam
│   │   ├── aln_trimmomatic_strictly_trimmed_sorted.bam.bai
│   │   ├── aln_trimmomatic_trimmed.bam
│   │   ├── aln_trimmomatic_trimmed.sam
│   │   ├── aln_trimmomatic_trimmed_sorted.bam
│   │   └── aln_trimmomatic_trimmed_sorted.bam.bai
│   ├── fastp_trimmed
│   │   ├── amp_res_R1_trimmed.fastq.gz
│   │   ├── amp_res_R2_trimmed.fastq.gz
│   │   ├── fastp_amp_res.html
│   │   └── fastp_amp_res.json
│   ├── fastqc
│   │   ├── amp_res_R1_fastqc.html
│   │   ├── amp_res_R1_fastqc.zip
│   │   ├── amp_res_R2_fastqc.html
│   │   └── amp_res_R2_fastqc.zip
│   ├── fastqc_fastp_trimmed
│   │   ├── amp_res_R1_trimmed_fastqc.html
│   │   ├── amp_res_R1_trimmed_fastqc.zip
│   │   ├── amp_res_R2_trimmed_fastqc.html
│   │   └── amp_res_R2_trimmed_fastqc.zip
│   ├── fastqc_trimmomatic_strictly_trimmed
│   │   ├── amp_res_strictly_1P_fastqc.html
│   │   ├── amp_res_strictly_1P_fastqc.zip
│   │   ├── amp_res_strictly_2P_fastqc.html
│   │   └── amp_res_strictly_2P_fastqc.zip
│   ├── fastqc_trimmomatic_trimmed
│   │   ├── amp_res_1P_fastqc.html
│   │   ├── amp_res_1P_fastqc.zip
│   │   ├── amp_res_2P_fastqc.html
│   │   └── amp_res_2P_fastqc.zip
│   ├── snpEff
│   │   ├── data
│   │   │   └── k12
│   │   │       ├── genes.gbk
│   │   │       ├── sequence.NC_000913.3.bin
│   │   │       └── snpEffectPredictor.bin
│   │   └── snpEff.config
│   ├── trimmomatic_strictly_trimmed
│   │   ├── amp_res_strictly_1P.fastq.gz
│   │   ├── amp_res_strictly_1U.fastq.gz
│   │   ├── amp_res_strictly_2P.fastq.gz
│   │   └── amp_res_strictly_2U.fastq.gz
│   ├── trimmomatic_trimmed
│   │   ├── amp_res_1P.fastq.gz
│   │   ├── amp_res_1U.fastq.gz
│   │   ├── amp_res_2P.fastq.gz
│   │   └── amp_res_2U.fastq.gz
│   └── variant_calling
│       ├── aln_fastp.mpileup
│       ├── aln_trimmomatic_strictly_trimmed.mpileup
│       ├── aln_trimmomatic_trimmed.mpileup
│       ├── varscan_res_aln_fastp.vcf
│       ├── varscan_res_aln_fastp_annotated.vcf
│       ├── varscan_res_aln_trimmomatic_strictly_trimmed.vcf
│       ├── varscan_res_aln_trimmomatic_strictly_trimmed_annotated.vcf
│       ├── varscan_res_aln_trimmomatic_trimmed.vcf
│       └── varscan_res_aln_trimmomatic_trimmed_annotated.vcf
├── scripts
│   ├── 01_downloading_raw_data_md5_check.sh
│   └── 02_collecting_aln_stats.sh
├── snpEff_genes.txt
└── snpEff_summary.html

18 directories, 100 files



################################################################
################ 13/11/2025 ####################################
################################################################

I tried to install MultiQC to bioinf environment, but that didn't work out.
MultiQC requires Python 3.11 and it doesn't support my curren Python 3.13.

```bash
mamba install -c bioconda multiqc 
```
(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ mamba install -c bioconda multiqc
warning  libmamba 'repo.anaconda.com', a commercial channel hosted by Anaconda.com, is used.
...
    └─ pin on python =3.13 * is not installable because it requires
       └─ python =3.13 *, which conflicts with any installable versions previously reported.
critical libmamba Could not solve for environment specs

To shortcut it, I installed MultiQC to its own mamba environment.
```bash
mamba create -n multiqc -c conda-forge -c bioconda python=3.11 multiqc
```

(bioinf) ➜  01_prac_ampicillin_resistance git:(master) ✗ mamba activate multiqc
(multiqc) ➜  01_prac_ampicillin_resistance git:(master) ✗ multiqc --version
multiqc, version 1.32

For some reason MultiQC don't see part of my fastqc reports. I dropped work with MultiQC for now.