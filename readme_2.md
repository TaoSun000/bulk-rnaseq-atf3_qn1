# 🧬 RNA-Seq Analysis: ATF3 Knockout vs Control SMCs

This repository contains the full workflow and outputs for a bulk RNA-seq analysis comparing **ATF3 knockout** and **control smooth muscle cells (SMCs)**, based on **NCBI BioProject [PRJNA716327](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA716327)**.

---

## 📁 Contents

- `RNA_Seq_analysis_question1_package.tar` — full source code, intermediate files, and final outputs
- `question_1_report.docx` — a formal report describing methods, results, and interpretations
- `README.md` — this documentation file

---

## 🔍 Objective

Compare gene expression between:
- **KO group**: SAMN18442667, SAMN18442665
- **Control group**: SAMN18442666, SAMN18442664  
to identify differentially expressed genes (DEGs) and enriched biological processes in **ATF3-deficient SMCs**.

---

## 🔧 Workflow Summary

| Step | Description | Tool Used |
|------|-------------|-----------|
| 1️⃣  | Download FASTQ from SRA | `prefetch`, `fasterq-dump` |
| 2️⃣  | Quality control | `FastQC` |
| 3️⃣  | Read trimming | `Trimmomatic` |
| 4️⃣  | Genome alignment | `HISAT2` (mm10 index) |
| 5️⃣  | Quantification | `featureCounts` |
| 6️⃣  | Differential expression | `DESeq2` in R |
| 7️⃣  | GO & KEGG enrichment | `clusterProfiler` in R |

---

## 📦 Software Versions

- **R**: 4.4.3  
- **DESeq2**: 1.40.2  
- **FastQC**: 0.11.9  
- **Trimmomatic**: 0.39  
- **HISAT2**: 2.1.0 (Ubuntu/WSL)  
- **Subread (featureCounts)**: 2.0.3  
- **clusterProfiler**: ≥4.8  
- **org.Mm.eg.db**: latest from Bioconductor

---

## 📊 Output Files

- `.sam`: alignment files for each sample
- `gene_counts.txt`: raw gene count matrix
- `DESeq2_results.csv`: log fold changes, p-values, and adjusted p-values
- `plots/`:  
  - MA plot  
  - Heatmap  
  - PCA plot  
  - GO barplot  
  - GO dotplot  
- `enrichment/`: output of GO and KEGG enrichment

---

## 🧪 Notable Findings

- **No DEGs** passed adjusted p-value < 0.05  
- Exploratory thresholding (unadjusted p < 0.05 & LFC > 1) yielded ~7 genes  
- GO enrichment revealed roles in:
  - calcium ion transport  
  - prostaglandin signaling  
  - epithelial morphogenesis  

---

## 📘 Citation

If you use this repository or pipeline, please cite:

- Love et al., *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*, Genome Biology (2014)
- Yu et al., *clusterProfiler: universal enrichment tool*, OMICS (2012)

---

## 🙋 Contact

Created by **Tao Sun**  
📧 Email: taosun618@gmail.com  
📍Ontario, Canada

---

