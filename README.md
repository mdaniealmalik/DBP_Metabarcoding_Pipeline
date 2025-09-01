# DBP_pipeline_nanopore_metabarcoding

# Introduction

The **DBP Metabarcoding Pipeline v2** is designed for processing Nanopore eDNA amplicon sequencing data.  
It performs the following steps:

1. **Read filtering** (NanoFilt)  
2. **Primer trimming** (Cutadapt)  
3. **OTU clustering** (VSEARCH)  
4. **Taxonomic assignment** (BLASTn)  
5. **Summary reporting**  

---

# Installation

Clone the repository and install dependencies using Conda:

```bash
conda env create -f DBP_dependencies_pipeline_v2.yml
conda activate DBP_pipeline
```
