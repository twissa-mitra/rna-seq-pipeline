# RNA-seq Mini-Pipeline (QC → Alignment → Quantification)

This repo demonstrates a minimal, reproducible RNA-seq pipeline using **WSL/Ubuntu** and standard tools.

It covers:

- **QC** with FastQC  
- **Alignment** with HISAT2  
- **Quantification** with featureCounts  
- **Basic sanity checks** with samtools

**Status:** pipeline runs end-to-end on a toy yeast dataset.  
**Next:** swap in a confirmed yeast RNA-seq dataset from GEO/SRA to get stronger counts.

---

## Quick Run (demo workflow)

All commands are run inside Ubuntu (WSL):

```bash
# enter project
cd ~/rnaseq

# align
hisat2 -x reference/sacCer3_index \
  -U data/yeast_subset.fastq \
  -S results/yeast_aligned.sam \
  --summary-file results/yeast_aligned.align.txt

# convert + sort + index
samtools view -bS results/yeast_aligned.sam | samtools sort -o results/yeast_aligned.sorted.bam
samtools index results/yeast_aligned.sorted.bam

# quantify
featureCounts -T 4 -a reference/sacCer3.gtf.gz -t exon -g gene_id \
  -M -O --fraction \
  -o results/gene_counts_exon_lenient.txt \
  results/yeast_aligned.sorted.bam


