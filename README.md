# DBP_pipeline_nanopore_metabarcoding

# Introduction

The **DBP Metabarcoding Pipeline v2** is designed for processing Nanopore eDNA amplicon sequencing data.  
It performs the following steps:

1. **Read filtering** (NanoFilt)  
2. **Primer trimming** (Cutadapt)  
3. **OTU clustering** (VSEARCH)  
4. **Taxonomic assignment** (BLASTn)  

---

# Installation

Clone the repository and install dependencies using Conda:

```bash
conda env create -f DBP_dependencies_pipeline_v2.yml
conda activate DBP_pipeline
```
Make the pipeline script executable:
```bash
chmod +x DBP_run_pipeline_v2.sh
```

📁 Project Folder Structure

When preparing your files for the pipeline, your working directory should be organized as follows:
MyProject/

```graphql
├── 1_Sample/                 # Raw input FASTQ files (gzip-compressed)
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── sample3.fastq.gz
│
├── database/                 # Reference database for BLAST
│   └── database.fasta        # Custom or curated reference sequences
│
├── DBP_run_pipeline_v2        # Pipeline script (make sure it's executable)
│
├── 2_NanoFilt_output/        # (Auto-generated) Quality filtered reads
├── 3_cutadapt_output/        # (Auto-generated) Primer-trimmed reads
├── 4_combined_fasta/         # (Auto-generated) FASTA converted reads
├── 5_vsearch/                # (Auto-generated) VSEARCH results
│   ├── rename_fasta/
│   └── combine/
│
├── otu_table.tsv             # (Auto-generated) OTU abundance table
└── result_blastn.txt         # (Auto-generated) BLAST taxonomic assignments
```

# Usage

Run with default parameters:
```bash
bash DBP_run_pipeline_v2.sh
```

Or customize parameters:
```bash
bash DBP_run_pipeline_v2.sh \
  -q 12 -l 200 -L 320 \
  --primer-fwd "TTTCTGTTGGTGCTGATATTGCGCCGGTAAAACTCGTGCCAGC" \
  --primer-rev "ACTTGCCTGTCGCTCTATCTTCCATAGTGGGGTATCTAATCCCAGTTTG" \
  --cutadapt-error 0.2 \
  --cutadapt-minlen 160 \
  --cutadapt-maxlen 220 \
  --vsearch-id 0.97 \
  --blast-evalue 1e-5 \
  --blast-identity 93 \
  --blast-qcov 95 \
  -t 8
```
PS: Sometimes the native barcode includes an overhang with the primer, so make sure to input the primer sequence with the overhang.

Parameter information:
| Tool         | Parameter           | Description                    |
| ------------ | ------------------- | ------------------------------ |
| **NanoFilt** | `-q 12`             | Min Phred quality score        |
|              | `-l 200`            | Min read length                |
|              | `-L 320`            | Max read length                |
| **Cutadapt** | `--primer-fwd`      | Forward primer sequence        |
|              | `--primer-rev`      | Reverse primer sequence        |
|              | `--cutadapt-error`  | Max allowed primer mismatch    |
|              | `--cutadapt-minlen` | Min read length after trimming |
|              | `--cutadapt-maxlen` | Max read length after trimming |
| **VSEARCH**  | `--vsearch-id`      | OTU clustering threshold       |
| **BLASTn**   | `--blast-evalue`    | E-value cutoff for hits        |
|              | `--blast-identity`  | Min percent identity           |
|              | `--blast-qcov`      | Min query coverage             |
| **CPU Threads**  | `-t`                | Number of CPU threads to use   |
