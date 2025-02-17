#!/bin/bash

# 📌 Définition des paramètres
BFILE="/home/p0129674/Documents/Emile_Analysis/Data/TZ.HbF.phased.3"
OUTDIR="/home/p0129674/Documents/Emile_Analysis/Data/resultats/qc_results"
PHENO_FILE="/home/p0129674/Documents/Emile_Analysis/regenie/pheno.txt"
#COVAR_FILE="covariates.txt"
BATCH_SIZE=1000
THREADS=8 # Nombre de processus parallèles

# 🧬 Liste des chromosomes
CHROMS=($(seq 1 22))

# 📂 Création du dossier de sortie
mkdir -p "$OUTDIR"

# 📢 Affichage des fichiers et paramètres
echo "📂 Fichier BED utilisé : $BFILE"
echo "🧬 Chromosomes analysés : ${CHROMS[@]}"
echo "📂 Résultats enregistrés dans : $OUTDIR"

# 🚀 Étape 1 : Contrôle de qualité avec PLINK
echo "🔍 Étape 1 : Contrôle de qualité avec PLINK..."
plink --bfile "$BFILE" --geno 0.02 --mind 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out "$OUTDIR/step1_filter"

plink --bfile "$OUTDIR/step1_filter" --indep-pairwise 50 5 0.2 --out "$OUTDIR/prune"
plink --bfile "$OUTDIR/step1_filter" --extract "$OUTDIR/prune.prune.in" --make-bed --out "$OUTDIR/step2_ld_pruned"

plink --bfile "$OUTDIR/step2_ld_pruned" --het --out "$OUTDIR/het_check"
awk '$6 > 0.05 || $6 < -0.05 {print $1, $2}' "$OUTDIR/het_check.het" > "$OUTDIR/remove_het_samples.txt"
plink --bfile "$OUTDIR/step2_ld_pruned" --remove "$OUTDIR/remove_het_samples.txt" --make-bed --out "$OUTDIR/QC_passed"

# 🚀 Étape 2 : Exécution de Regenie en parallèle pour chaque chromosome
echo "🚀 Étape 2 : Exécution de Regenie sur chaque chromosome en parallèle..."

run_regenie() {
    chrom=$1
    echo "🧬 Traitement du chromosome $chrom..."
    
    docker run --rm -v "$OUTDIR:/data" regenie regenie --step 1 \
        --bed /data/QC_passed \
        --phenoFile /data/$PHENO_FILE \
        --covarFile /data/$COVAR_FILE \
        --bsize $BATCH_SIZE \
        --out /data/regenie_step1_chr$chrom

    docker run --rm -v "$OUTDIR:/data" regenie regenie --step 2 \
        --bed /data/QC_passed \
        --phenoFile /data/$PHENO_FILE \
        --covarFile /data/$COVAR_FILE \
        --bsize $BATCH_SIZE \
        --pred /data/regenie_step1_chr${chrom}_pred.list \
        --chr $chrom \
        --out /data/regenie_chr$chrom.out
}

export -f run_regenie
parallel -j $THREADS run_regenie ::: "${CHROMS[@]}"

echo "✅ Analyse GWAS terminée !"
