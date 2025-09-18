#!/usr/bin/env bash
set -euo pipefail

# ---------- Config (override via env: THREADS, FASTQ, etc.) ----------
THREADS="${THREADS:-4}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${DATA_DIR:-$ROOT/data}"
REF_DIR="${REF_DIR:-$ROOT/ref}"
ALN_DIR="${ALN_DIR:-$ROOT/aln}"
QC_DIR="${QC_DIR:-$ROOT/qc}"
COUNTS_DIR="${COUNTS_DIR:-$ROOT/counts}"

FASTQ="${FASTQ:-$DATA_DIR/yeast2.fastq.gz}"
GENOME_FA="$REF_DIR/sacCer3.fa"
INDEX_PREFIX="$REF_DIR/sacCer3"
GTF="$REF_DIR/sacCer3.ensGene.gtf"
SGD_MAP_GZ="$REF_DIR/sgdToName.txt.gz"

# Yeast FASTQ (good SRR we used)
FASTQ_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR295/061/SRR29550561/SRR29550561.fastq.gz"
# UCSC sacCer3 reference + genes + SGD name map
FA_URL="https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz"
GTF_URL="https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz"
SGD_MAP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/sgdToName.txt.gz"

# ---------- Helpers ----------
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing tool: $1"; exit 127; }; }
for t in fastqc hisat2 hisat2-build samtools featureCounts multiqc awk zcat wget gzip; do need "$t"; done

mkdir -p "$DATA_DIR" "$REF_DIR" "$ALN_DIR" "$QC_DIR/multiqc" "$COUNTS_DIR"

# ---------- 1) Data ----------
if [[ ! -f "$FASTQ" ]]; then
  echo "[DATA] Downloading FASTQ …"
  wget -c -O "$FASTQ" "$FASTQ_URL"
fi
echo "[DATA] Verifying FASTQ gzip…"
gzip -t "$FASTQ" || { echo "FASTQ appears corrupt"; exit 1; }

echo "[QC] FastQC…"
fastqc "$FASTQ" -o "$QC_DIR"

# ---------- 2) Reference ----------
if [[ ! -f "$GENOME_FA" ]]; then
  echo "[REF] Downloading sacCer3 FASTA…"
  wget -c -O "$REF_DIR/sacCer3.fa.gz" "$FA_URL"
  gzip -t "$REF_DIR/sacCer3.fa.gz"
  zcat "$REF_DIR/sacCer3.fa.gz" > "$GENOME_FA"
fi

if [[ ! -e "${INDEX_PREFIX}.1.ht2" ]]; then
  echo "[REF] Building HISAT2 index…"
  hisat2-build -p "$THREADS" "$GENOME_FA" "$INDEX_PREFIX"
fi

if [[ ! -f "$GTF" ]]; then
  echo "[REF] Downloading GTF…"
  wget -c -O "$REF_DIR/sacCer3.ensGene.gtf.gz" "$GTF_URL"
  gzip -t "$REF_DIR/sacCer3.ensGene.gtf.gz"
  zcat "$REF_DIR/sacCer3.ensGene.gtf.gz" > "$GTF"
fi

if [[ ! -f "$SGD_MAP_GZ" ]]; then
  echo "[REF] Downloading SGD name map…"
  wget -c -O "$SGD_MAP_GZ" "$SGD_MAP_URL"
fi

# ---------- 3) Align ----------
SAM="$ALN_DIR/yeast2.sam"
BAM="$ALN_DIR/yeast2.sorted.bam"

echo "[ALIGN] HISAT2…"
hisat2 -p "$THREADS" -x "$INDEX_PREFIX" -U "$FASTQ" -S "$SAM" \
  --summary-file "$ALN_DIR/yeast2.hisat2.summary.txt"

echo "[ALIGN] Sort & index BAM…"
samtools view -@ "$THREADS" -bS "$SAM" | samtools sort -@ "$THREADS" -o "$BAM" -
samtools index "$BAM"
samtools flagstat "$BAM" | tee "$ALN_DIR/yeast2.flagstat.txt"
rm -f "$SAM"

# ---------- 4) Count ----------
echo "[COUNT] featureCounts…"
featureCounts -T "$THREADS" -a "$GTF" -t exon -g gene_id \
  -o "$COUNTS_DIR/yeast2.featureCounts.txt" "$BAM"

echo "[COUNT] Tidy counts + gene names…"
awk 'NR>2{print $1"\t"$7}' "$COUNTS_DIR/yeast2.featureCounts.txt" > "$COUNTS_DIR/yeast2.counts.tsv"
zcat "$SGD_MAP_GZ" | \
  awk 'NR==FNR{m[$1]=$2; next} {nm=( ($1 in m)?m[$1]:$1 ); print $1"\t"nm"\t"$2}' \
  - "$COUNTS_DIR/yeast2.counts.tsv" > "$COUNTS_DIR/yeast2.counts.with_names.tsv"
awk 'BEGIN{print "gene_id,gene_name,count"} {print $1","$2","$3}' \
  "$COUNTS_DIR/yeast2.counts.with_names.tsv" > "$COUNTS_DIR/yeast2.counts.with_names.csv"

# ---------- 5) MultiQC ----------
echo "[QC] MultiQC…"
multiqc "$ROOT" -o "$QC_DIR/multiqc"

echo "[DONE] Results:"
echo "  - MultiQC: $QC_DIR/multiqc/multiqc_report.html"
echo "  - Counts:  $COUNTS_DIR/yeast2.counts.with_names.csv"
