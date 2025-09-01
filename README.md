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

# Usage

Run with default parameters:
```bash
bash DBP_run_pipeline_v2.sh
```

Or customize parameters:
```bash
bash DBP_run_pipeline_v2.sh \
  -q 12 -l 200 -L 320 \
  --primer-fwd "GACTGGGATTAGATACCCC" \
  --primer-rev "TGGCACGAGTTTTACCGGC" \
  --cutadapt-error 0.1 \
  --cutadapt-minlen 160 \
  --cutadapt-maxlen 220 \
  --vsearch-id 0.97 \
  --blast-evalue 1e-5 \
  --blast-identity 93 \
  --blast-qcov 95 \
  -t 8
```

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
