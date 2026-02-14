# circRNA Detection and Differential Expression using CIRI3

This repository provides a complete and reproducible pipeline for **circRNA detection and differential expression analysis** using CIRI3.

The workflow includes:

* RNA-seq data download from SRA
* Genome indexing and alignment using BWA
* circRNA detection (single and multi-sample modes)
* Gene quantification using featureCounts
* Differential expression analysis:

  * **DE_BSJ** (absolute circRNA abundance)
  * **DE_Ratio** (circularization efficiency)
  * **DE_Relative** (isoform switching)

This pipeline is designed for paired-end RNA-seq data and is tested on dataset **GSE97239** (human, GRCh38).

---

## 📌 Overview of Workflow

```
SRA → FASTQ → Alignment (BWA-MEM) → CIRI3 detection
                              ↓
                    featureCounts (gene counts)
                              ↓
                Differential Expression (3 models)
```

---

## ⚙️ Environment Setup

### 1️⃣ Install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda --version
```

---

### 2️⃣ Clone CIRI3

```bash
mkdir -p ~/bioinfo_tools
cd ~/bioinfo_tools
git clone https://github.com/gyjames/CIRI3.git
cd CIRI3
```

---

### 3️⃣ Create Conda Environment

```bash
conda env create -n CIRI3 -f ./environment.yaml
conda activate CIRI3
```

Install required tools:

```bash
conda install -c bioconda bwa samtools sra-tools subread -y
conda install -c conda-forge openjdk=8 -y
```

Test installation:

```bash
java -version
bwa
samtools --version
featureCounts -v
```

---

## 📂 Project Structure

```
CIRI3_analysis/
├── data/
│   ├── raw_fastq/
│   └── reference/
├── results/
│   ├── alignment/
│   ├── ciri3/
│   └── DE/
└── logs/
```

---

## 📥 Data Preparation

### Reference Genome (GRCh38)

```bash
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip *.gz
```

Index genome:

```bash
bwa index -a bwtsw Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

---

### Download RNA-seq Data

Recommended (faster than fastq-dump):

```bash
fasterq-dump SRRXXXXXXX -S -e 8
gzip SRRXXXXXXX_*.fastq
```

---

## 🧬 Alignment (BWA-MEM)

Important: use `-T 19` to retain split alignments required for circRNA detection.

```bash
bwa mem -T 19 -t 8 \
    reference.fa \
    sample_1.fastq.gz sample_2.fastq.gz | \
    samtools sort -@ 8 -o sample.sorted.bam
```

---

## 🔎 circRNA Detection with CIRI3

### Single-Sample Mode (Required for DE_BSJ)

```bash
java -jar CIRI3_Java_1.8.0.jar \
    -I sample.sam \
    -O sample_ciri3 \
    -F genome.fa \
    -A annotation.gtf
```

---

### Multi-Sample Mode (Required for DE_Ratio and DE_Relative)

Create `sample_list.txt` with absolute paths to SAM files:

```bash
java -jar CIRI3_Java_1.8.0.jar \
    -I sample_list.txt \
    -O all_samples.txt \
    -F genome.fa \
    -A annotation.gtf \
    -W 1 -T 8
```

Generated outputs:

* `all_samples.txt`
* `all_samples.txt.BSJ_Matrix`
* `all_samples.txt.FSJ_Matrix`

---

## 🧮 Gene Expression Quantification

```bash
featureCounts -p -T 8 -t exon -g gene_id \
    -a annotation.gtf \
    -o gene_counts.txt \
    *.sorted.bam
```

Format gene expression matrix for CIRI3 compatibility before DE_BSJ.

---

# 📊 Differential Expression Models

CIRI3 provides three complementary DE methods:

---

## 1️⃣ DE_BSJ – Absolute circRNA Abundance

Tests overall circRNA expression changes.

```bash
java -jar CIRI3_Java_1.8.0.jar DE_BSJ \
    -I sample_info_bsj.tsv \
    -G gene_expression.txt \
    -O BSJ_DE_results.txt \
    -P 0.05
```

✔ Detects circRNAs up/downregulated between conditions.

---

## 2️⃣ DE_Ratio – Circularization Efficiency

Tests:

```
BSJ / (BSJ + FSJ)
```

```bash
java -jar CIRI3_Java_1.8.0.jar DE_Ratio \
    -I sample_info_ratio.tsv \
    -BM all_samples.txt.BSJ_Matrix \
    -FM all_samples.txt.FSJ_Matrix \
    -O JR_DE_results.txt \
    -T 8
```

✔ Identifies changes in circular-to-linear transcript ratio.

---

## 3️⃣ DE_Relative – Isoform Switching

Detects shifts in relative abundance of circRNA isoforms from the same host gene.

```bash
java -jar CIRI3_Java_1.8.0.jar DE_Relative \
    -I sample_info_ratio.tsv \
    -M all_samples.txt.BSJ_Matrix \
    -GC circ_gene.txt \
    -O RE_switching_results.txt \
    -T 8
```

✔ Reveals circRNA isoform switching events.

---

# 🧠 Biological Interpretation

| Model       | Tests                      | Interpretation              |
| ----------- | -------------------------- | --------------------------- |
| DE_BSJ      | Absolute abundance         | Overall circRNA regulation  |
| DE_Ratio    | Circularization efficiency | Splicing regulation changes |
| DE_Relative | Isoform proportion shifts  | Isoform switching           |

Using all three models provides a systems-level understanding of circRNA regulation.

---

# 🧪 Recommended Filtering

Before interpretation:

* BSJ count ≥ 2–5 in at least 2 samples
* Remove low-expression circRNAs
* Adjust for multiple testing (FDR)

---

# 📦 Software Versions

* CIRI3 (Java 1.8)
* BWA-MEM
* SAMtools
* featureCounts (Subread)
* edgeR (via CIRI3)
