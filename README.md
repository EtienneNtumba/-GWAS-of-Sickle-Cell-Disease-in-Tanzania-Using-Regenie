# **Genome-Wide Association Study (GWAS) of Sickle Cell Disease in Tanzania Using Regenie**

## **1. Introduction**
Genome-wide association studies (GWAS) are used to identify genetic variants associated with complex traits and diseases. This study focuses on **Sickle Cell Disease (SCD) in Tanzania**, aiming to detect **genetic variants significantly associated with the phenotype of interest** using **Regenie**, a two-step approach suitable for large-scale GWAS.

**Key Details:**
- **Total Samples:** 3,210 individuals.
- **Original SNPs:** 8,457,145 (before quality control).
- **Post-QC SNPs:** 1,466,733 (after applying quality control).
- **Phenotype Type:** Continuous variable.
- **Software Used:** `PLINK` for quality control, `Regenie` for association testing.
- **Analysis Type:** GWAS with **stepwise regression modeling in Regenie**.

---

## **2. Quality Control (QC) Using PLINK**
Before conducting a GWAS, it is **crucial to perform QC** to ensure **reliable** and **valid** results by removing:
- **Low-quality SNPs** (e.g., missing data, low minor allele frequency).
- **Samples with excess heterozygosity** (potential genotyping errors).
- **Population structure issues** (through pruning of correlated SNPs).

### **Step 1: SNP and Sample Filtering**
We apply the following **PLINK filters**:

```bash
plink --bfile "$BFILE" \
      --geno 0.02 \  # Remove SNPs with >2% missing genotypes
      --mind 0.02 \  # Remove individuals with >2% missing genotypes
      --maf 0.01 \   # Exclude SNPs with Minor Allele Frequency (MAF) <1%
      --hwe 1e-6 \   # Hardy-Weinberg Equilibrium filter (p < 1e-6)
      --make-bed --out "$OUTDIR/step1_filter"
```

**Why?**
- **`--geno 0.02`** ensures only SNPs with low missingness are kept.
- **`--mind 0.02`** removes individuals with excessive missing genotypes.
- **`--maf 0.01`** keeps SNPs with sufficient variability for analysis.
- **`--hwe 1e-6`** removes SNPs deviating from Hardy-Weinberg equilibrium in controls (suggests possible genotyping errors).

### **Step 2: Pruning for Linkage Disequilibrium (LD)**
To avoid **correlated SNPs affecting association results**, we **prune SNPs** based on LD.

```bash
plink --bfile "$OUTDIR/step1_filter" \
      --indep-pairwise 50 5 0.2 \  # Window size 50, step 5, R^2 threshold 0.2
      --out "$OUTDIR/prune"

plink --bfile "$OUTDIR/step1_filter" \
      --extract "$OUTDIR/prune.prune.in" \
      --make-bed --out "$OUTDIR/step2_ld_pruned"
```

**Why?**
- LD pruning **removes redundant SNPs** that are highly correlated.
- Keeps a **set of independent SNPs**, reducing redundancy.

### **Step 3: Heterozygosity Check**
Heterozygosity outliers indicate **sample contamination or genotyping errors**.

```bash
plink --bfile "$OUTDIR/step2_ld_pruned" --het --out "$OUTDIR/het_check"

awk '$6 > 0.05 || $6 < -0.05 {print $1, $2}' "$OUTDIR/het_check.het" > "$OUTDIR/remove_het_samples.txt"

plink --bfile "$OUTDIR/step2_ld_pruned" \
      --remove "$OUTDIR/remove_het_samples.txt" \
      --make-bed --out "$OUTDIR/QC_passed"
```

**Why?**
- Samples with **excess heterozygosity (|F coefficient| > 0.05)** may be **misgenotyped or contaminated**.

---

## **3. GWAS Analysis Using Regenie**
### **Why Regenie?**
- Efficient for **large-scale genetic data**.
- Uses **stepwise ridge regression**, reducing confounding by population structure.

  # 🧬 GWAS Analysis with REGENIE - Installation & Visualization Guide

This repository provides a **step-by-step guide** on how to **install, run, and visualize GWAS results using REGENIE**.

---

## 🚀 3.1. Installing REGENIE

You can install **REGENIE** using **Docker** (recommended) or **Conda**.

### 📌 Option 1: Install REGENIE using Docker
Docker is the easiest method as it avoids dependency issues.

#### 🔹 Step 1: Pull the REGENIE Docker Image
```bash
docker pull rgcgithub/regenie

docker run --rm rgcgithub/regenie regenie --help

## 📌 Option 2: Install REGENIE using Conda

If you prefer to install REGENIE locally, use Conda.

🔹 Step 1: Create a Conda Environment

conda create -n regenie_env -y
conda activate regenie_env

🔹 Step 2: Install Dependencies

conda install -c conda-forge cmake make gcc libgomp htslib boost -y

🔹 Step 3: Download and Compile REGENIE

git clone https://github.com/rgcgithub/regenie.git
cd regenie
make

🔹 Step 4: Verify the Installation

./regenie --help


### **Step 1: Null Model Fitting**
```bash
docker run --rm -v ~/Documents/Emile_Analysis/Data/resultats/qc_results:/data -w /data regenie regenie \
  --step 1 \
  --bed QC_passed_filtered \
  --phenoFile pheno.txt \
  --covarFile covariates.txt \
  --lowmem \
  --bsize 1000 \
  --out step1
```

**Key Parameters:**
- **`--step 1`**: Builds the null model (ridge regression).
- **`--lowmem`**: Optimizes memory for large datasets.
- **`--bsize 1000`**: Block size for model fitting.

### **Step 2: Association Testing**
```bash
docker run --rm -v ~/Documents/Emile_Analysis/Data/resultats/qc_results:/data -w /data regenie regenie \
  --step 2 \
  --bed QC_passed \
  --phenoFile pheno.txt \
  --covarFile covariates.txt \
  --pred step1_pred.list \
  --bsize 1000 \
  --out step2
```

---

## **4. Results & Interpretation**
### **4.1 Top SNPs (p-value < 1e-4)**
The top SNPs from the analysis indicate significant associations.

| CHROM | GENPOS     | BETA   | SE     | LOG10P  |
|-------|-----------|--------|--------|---------|
| 1     | 6066007   | 5.4914 | 1.2504 | 4.9489  |
| 1     | 7284761   | 7.0988 | 1.7273 | 4.4024  |
| 1     | 21731734  | 4.5567 | 1.1240 | 4.2984  |
| 1     | 52065227  | 9.2386 | 2.2894 | 4.2635  |
| 1     | 208705845 | 7.6419 | 1.8555 | 4.4188  |
| 1     | 208775262 | 9.4583 | 2.3513 | 4.2400  |
| 1     | 244358190 | 5.2966 | 1.3178 | 4.2341  |
| 2     | 77757146  | 6.0016 | 1.4618 | 4.3945  |
| 2     | 77772143  | 6.0016 | 1.4618 | 4.3945  |
| 2     | 164577126 | 10.1276| 2.5565 | 4.1280  |
| 2     | 181776641 | 9.5743 | 2.4598 | 4.0031  |
| 2     | 185580674 | 7.6699 | 1.8963 | 4.2808  |
| 2     | 199114681 | 6.8673 | 1.6646 | 4.4318  |
| 2     | 216545159 | 5.9755 | 1.5010 | 4.1634  |
| 3     | 2370966   | 7.8825 | 1.9740 | 4.1858  |
| 3     | 5391627   | 5.6195 | 1.3944 | 4.2533  |
| 3     | 60304669  | 5.5119 | 1.4164 | 4.0016  |
| 3     | 60403966  | 6.2029 | 1.4050 | 4.9952  |
| 3     | 60412826  | 5.9351 | 1.5083 | 4.0799  |
| 3     | 60432259  | 5.9351 | 1.5083 | 4.0799  |
| 3     | 60464649  | 7.8818 | 1.9739 | 4.1856  |
| 3     | 73377487  | 5.0682 | 1.2899 | 4.0694  |
| 3     | 98121155  | 8.6099 | 2.1980 | 4.0476  |
| 3     | 100240260 | 7.1867 | 1.7437 | 4.4246  |
| 3     | 154899364 | 9.5532 | 2.3445 | 4.3366  |
| 4     | 31255542  | 4.2436 | 1.0753 | 4.1005  |
| 11    | 79615053  | 3.7527 | 0.9586 | 4.0433  |
| 13    | 74663941  | 5.4198 | 1.3156 | 4.4209  |
| 14    | 66420147  | 9.0633 | 2.2660 | 4.1977  |
| 14    | 82413564  | 5.2212 | 1.2803 | 4.3427  |
| 14    | 89985346  | 7.1903 | 1.7660 | 4.3305  |
| 14    | 100566037 | 8.2298 | 2.0887 | 4.0891  |
| 14    | 101533093 | 8.9655 | 2.2725 | 4.0983  |
| 14    | 105373185 | 7.6832 | 1.9133 | 4.2272  |
| 15    | 36704126  | 4.3857 | 1.0854 | 4.2731  |
| 18    | 35573548  | 5.5180 | 1.3435 | 4.3973  |
| 18    | 77361904  | 4.7720 | 1.2124 | 4.0817  |



### **4.2 QQ Plot**
![QQ Plot](qq_plot.png)

### **4.3 Manhattan Plot**
![Manhattan Plot](manhattan_plot.png)

---

## **5. Next Steps & Improvements**
- **Validate Significant SNPs** using fine-mapping.
- **Replication in Independent Cohorts**.
- **Functional Annotation** using Annovar or RegulomeDB.
- **Alternative Methods** like SAIGE or BOLT-LMM.

---

## **6. Conclusion**
This GWAS study successfully identified **significant SNPs** for **Sickle Cell Disease in Tanzania**, demonstrating the effectiveness of **Regenie** in handling large-scale genetic data.

---

