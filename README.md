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



# 💻 Installation (Step by Step)

Clone the repository and install dependencies using Conda:

**Clone the repository**  
```bash
git clone https://github.com/mdaniealmalik/DBP_pipeline_nanopore_metabarcoding.git
```
```bash
cd DBP_pipeline_nanopore_metabarcoding
```

**Create and activate the conda environment**
```bash
conda env create -f  environment.yml
```
```bash
conda activate dbp_pipeline
```

**Make the pipeline script executable**
```bash
chmod +x bin/DBP_run_pipeline_v2.sh
```

# 📁 Project Folder Structure

When preparing your files for the pipeline, your working directory should be organized as follows:

```graphql
MyProject/
├── 1_Sample/                 # Raw input FASTQ files (gzip-compressed)
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── sample3.fastq.gz
│
├── database/                 # Reference database for BLAST
│   └── database.fasta        # Custom or curated reference sequences (only accession number at the header)
│   └── database.txt          # Custom or curated reference taxon names (Accession number and taxon information) 
│
├── DBP_run_pipeline_v2.sh    # Pipeline script (make sure it's executable)
├── DBP_LCA_assign.R          # Pipeline script for LCA (make sure it's executable)
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

You can use the data structure in [my example data and database](https://github.com/mdaniealmalik/DBP_pipeline_nanopore_metabarcoding/tree/main/MyProject) to try this pipeline, and before trying it on your real datasets. 

## 📚 Database
This pipeline requires a reference database for taxonomic assignment with **BLASTn**.  
Place your reference FASTA file in the folder `database/`:  
- The pipeline will automatically create BLAST indices (`.nin`, `.nhr`, `.nsq`) from this FASTA file.  

You may use:
- A curated database (e.g., **MitoFish 12S**, **NCBI RefSeq**)  
- A custom database built from your target taxa

🐟 We have prepared a curated sequence database (Database Ikan Indonesia or DAKI) for a data sequence fit for MiFish, containing only marine fish species from Indonesia. Please visit [this link](https://github.com/mdaniealmalik/DBP_pipeline_nanopore_metabarcoding/tree/main/Curated-Metabarcoding-Database-for-Marine-Fish-in-Indonesia) to download the database and view details of the curation process.

⚠️ Make sure the file is named `database.fasta` and stored inside the `database/` folder before running the pipeline.

# ▶️ Usage
Before running the pipeline on your dataset, copy the file `DBP_run_pipeline_v2.sh` and `DBP_LCA_assign.R` from the `bin` folder into your dataset folder (e.g., `MyProject`), and then navigate into the `MyProject` directory.

```bash
cp bin/DBP_run_pipeline_v2.sh MyProject/DBP_run_pipeline_v2.sh
```
```bash
cp bin/DBP_LCA_assign.R MyProject/DBP_LCA_assign.R
```

```bash
cd MyProject
```
Run with default parameters:
```bash
bash DBP_run_pipeline_v2.sh
```

Or customise parameters:
```bash
bash DBP_run_pipeline_v2.sh \
  -q 12 -l 200 -L 320 \
  --primer-fwd "TTTCTGTTGGTGCTGATATTGCGCCGGTAAAACTCGTGCCAGC" \
  --primer-rev "ACTTGCCTGTCGCTCTATCTTCCATAGTGGGGTATCTAATCCCAGTTTG" \
  --cutadapt-error 0.2 \
  --cutadapt-minlen 160 \
  --cutadapt-maxlen 220 \
  --vsearch-id 0.97 \
  --blast-evalue 0.001 \
  --blast-identity 93 \
  --blast-qcov 95 \
  --blast_max_target 25\
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
|              | `--blast_max_target`          | how many top database sequences (hits) (**Default: 25**)  |
| **CPU Threads**  | `-t`                | Number of CPU threads to use (**Default: 4**)  |


**Run the optional LCA refinement**
```bash
Rscript DBP_LCA_assign.R
```
Or customise parameter (follow your minimum identity from blastn):
```bash
Rscript DBP_LCA_assign.R 97
```
