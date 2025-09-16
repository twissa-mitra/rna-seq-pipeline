# RNA-seq Pipeline (QC → Alignment → Quantification → DE)

This repo demonstrates a minimal, reproducible RNA-seq pipeline.

It covers:

\- QC with FastQC

\- Alignment with HISAT2

\- Quantification with featureCounts

\- Differential expression with DESeq2



\*\*Status:\*\* setting up environment \& data.  

\*\*Next:\*\* download a small public dataset (GEO/SRA), run QC and alignment, push results.





---



\## Toy DESeq2 Demo (quick check)



This repo includes a tiny, self-contained DESeq2 demo (no real data) to verify the R side works.



\*\*Run it:\*\*

1\. Open `notebooks/deseq2\_toy.R` in RStudio.

2\. In RStudio: `Session ▸ Set Working Directory ▸ To Source File Location`.

3\. Click \*\*Source\*\* to run the whole script.



\*\*Outputs go to `results/`:\*\*

\- `deseq2\_toy\_results.csv` — toy differential results  

\- `deseq2\_toy\_normalized\_counts.csv` — normalized toy counts  

\- `deseq2\_toy\_MAplot.pdf` — MA plot image  

