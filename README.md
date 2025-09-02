# DBP_pipeline_nanopore_metabarcoding

# 📖 Introduction

The **DBP Metabarcoding Pipeline v2** is designed for Nanopore eDNA amplicon sequencing data.  
Please note that before using this pipeline, it is assumed that Dorado (or another Nanopore basecaller) has already performed:

- ✅ Basecalling  
- ✅ Demultiplexing (per-sample FASTQ files)  
- ✅ Adapter and barcode removal  

---

## 🧬 Pipeline Workflow  

The pipeline then proceeds with the following steps:

1. **Quality Filtering** – remove low-quality or too short/long reads  
   - Tool: `NanoFilt`  

2. **Primer Trimming** – remove amplification primers  
   - Tool: `Cutadapt`  

3. **OTU Clustering** – dereplicate, cluster, and remove chimeras  
   - Tool: `VSEARCH`  

4. **Taxonomic Assignment** – assign sequences to reference database  
   - Tool: `BLASTn`  

5. **Results** – generate:  
   - `otu_table.tsv` (OTU abundance table)  
   - `result_blastn.txt` (BLAST taxonomy results)  



# 💻 Installation

Clone the repository and install dependencies using Conda:

```bash
conda env create -f DBP_dependencies_pipeline_v2.yml
conda activate DBP_pipeline
```
Make the pipeline script executable:
```bash
chmod +x DBP_run_pipeline_v2.sh
```

# 📁 Project Folder Structure

When preparing your files for the pipeline, your working directory should be organized as follows:
MyProject/

```graphql
├── 1_Sample/                 # Raw input FASTQ files (gzip-compressed)
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── sample3.fastq.gz
│
├── database/                 # Reference database for BLAST
│   └── database.fasta        # Custom or curated reference sequences (only accession number at the header)
│   └── database.txt          # Custom or curated reference taxon names (Accesion number and taxon information) 
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

# ▶️ Usage

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
**Note**: Sometimes the native barcode includes an overhang with the primer, so make sure to input the primer sequence with the overhang.

Parameter information:
| Tool         | Parameter           | Description                    |
| ------------ | ------------------- | ------------------------------ |
| **NanoFilt** | `-q`             | Min Phred quality score (**Default: 12**)       |
|              | `-l`            | Min read length (**Default: 180**)               |
|              | `-L`            | Max read length (**Default: 320**)               |
| **Cutadapt** | `--primer-fwd`      | Forward primer sequence  (**Default: overhang + Mifish-U**)      |
|              | `--primer-rev`      | Reverse primer sequence  (**Default: overhang + Mifish-U**)      |
|              | `--cutadapt-error`  | Max allowed primer mismatch  (**Default: 0.2**)  |
|              | `--cutadapt-minlen` | Min read length after trimming (**Default: 150**) |
|              | `--cutadapt-maxlen` | Max read length after trimming (**Default: 200**) |
| **VSEARCH**  | `--vsearch-id`      | OTU clustering threshold (**Default: 97**)      |
| **BLASTn**   | `--blast-evalue`    | E-value cutoff for hits (**Default: 0.001**)    |
|              | `--blast-identity`  | Min percent identity (**Default: 0.97**)         |
|              | `--blast-qcov`      | Min query coverage  (**Default: 0.90**)         |
| **CPU Threads**  | `-t`                | Number of CPU threads to use (**Default: 4**)  |
