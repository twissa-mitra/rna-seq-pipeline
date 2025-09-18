# Yeast RNA-seq Mini Pipeline (QC → Align → Count)

A tiny, reproducible RNA-seq example using *Saccharomyces cerevisiae* (sacCer3).

## What this shows
- Download & verify FASTQ
- FastQC quality check
- Adapter trim (if needed) with cutadapt
- Align with HISAT2 to *sacCer3*
- Count genes with featureCounts
- Summarize with MultiQC
- Deliver tidy counts with gene names

## Quick start (commands actually run)
```bash
# Data
wget -O data/yeast2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR295/061/SRR29550561/SRR29550561.fastq.gz
gunzip -t data/yeast2.fastq.gz && echo OK

# QC
fastqc data/yeast2.fastq.gz -o qc/

# Reference (UCSC sacCer3)
wget -O ref/sacCer3.fa.gz https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gzip -t ref/sacCer3.fa.gz && zcat ref/sacCer3.fa.gz > ref/sacCer3.fa
hisat2-build -p 4 ref/sacCer3.fa ref/sacCer3

# Align
hisat2 -p 4 -x ref/sacCer3 -U data/yeast2.fastq.gz -S aln/yeast2.sam --summary-file aln/yeast2.hisat2.summary.txt
samtools view -@ 4 -bS aln/yeast2.sam | samtools sort -@ 4 -o aln/yeast2.sorted.bam -
samtools index aln/yeast2.sorted.bam
samtools flagstat aln/yeast2.sorted.bam | tee aln/yeast2.flagstat.txt

# Annotation (UCSC Ensembl genes + SGD names)
wget -O ref/sacCer3.ensGene.gtf.gz https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
gzip -t ref/sacCer3.ensGene.gtf.gz && zcat ref/sacCer3.ensGene.gtf.gz > ref/sacCer3.ensGene.gtf
wget -O ref/sgdToName.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/sgdToName.txt.gz

# Counts
featureCounts -T 4 -a ref/sacCer3.ensGene.gtf -t exon -g gene_id \
  -o counts/yeast2.featureCounts.txt aln/yeast2.sorted.bam

# Tidy counts
awk 'NR>2{print $1"\t"$7}' counts/yeast2.featureCounts.txt > counts/yeast2.counts.tsv
zcat ref/sgdToName.txt.gz | \
awk 'NR==FNR{m[$1]=$2; next} {nm=( ($1 in m)?m[$1]:$1 ); print $1"\t"nm"\t"$2}' \
  - counts/yeast2.counts.tsv > counts/yeast2.counts.with_names.tsv
awk 'BEGIN{print "gene_id,gene_name,count"} {print $1","$2","$3}' \
  counts/yeast2.counts.with_names.tsv > counts/yeast2.counts.with_names.csv

# MultiQC summary
multiqc . -o qc/multiqc



## Install & Run

### 1) Install tools (Ubuntu/WSL)
```bash
bash scripts/install_tools.sh
```

### 2) Run the pipeline
```bash
THREADS=4 bash scripts/run_pipeline.sh
```

### Outputs
- `qc/multiqc/multiqc_report.html`
- `counts/yeast2.counts.with_names.csv`
- Tool versions: `env/versions.txt`
