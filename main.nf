#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 📂 Définition des paramètres
params.bfile   = "/home/p0129674/Documents/Emile_Analysis/Data/TZ.HbF.phased.3"
params.outdir  = "/home/p0129674/Documents/Emile_Analysis/Data/resultats/qc_results"
params.chroms  = (1..22) // 🧬 Liste des chromosomes

// 📢 Affichage des fichiers et chromosomes détectés
println "📂 Fichiers BED trouvés: ${params.bfile}"
println "🧬 Chromosomes traités: ${params.chroms}"

// 🔹 Création des channels
channel_chroms = channel.fromList(params.chroms)
channel_files = channel.from([
    file("${params.bfile}.bed"),
    file("${params.bfile}.bim"),
    file("${params.bfile}.fam")
])

// 🔹 Processus de Contrôle de Qualité avec PLINK
process quality_control {
    container 'plink'

    input:
    path bed
    path bim
    path fam

    output:
    path "QC_passed.bed", emit: bed
    path "QC_passed.bim", emit: bim
    path "QC_passed.fam", emit: fam

    script:
    """
    plink --bfile ${bed.baseName} --geno 0.02 --mind 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out step1_filter
    plink --bfile step1_filter --indep-pairwise 50 5 0.2 --out prune
    plink --bfile step1_filter --extract prune.prune.in --make-bed --out step2_ld_pruned
    plink --bfile step2_ld_pruned --het --out het_check
    awk '\$6 > 0.05 || \$6 < -0.05 {print \$1, \$2}' het_check.het > remove_het_samples.txt
    plink --bfile step2_ld_pruned --remove remove_het_samples.txt --make-bed --out QC_passed
    """
}

// 🔹 Processus d'exécution de Regenie (par chromosome)
process run_regenie {
    container 'regenie'

    input:
    path bed
    path bim
    path fam
    val chrom from channel_chroms

    output:
    path "regenie_chr${chrom}.out"

    script:
    """
    regenie --step 1 \
        --bed QC_passed \
        --phenoFile pheno.txt \
        --covarFile covariates.txt \
        --bsize 1000 \
        --out regenie_step1_chr${chrom}

    regenie --step 2 \
        --bed QC_passed \
        --phenoFile pheno.txt \
        --covarFile covariates.txt \
        --bsize 1000 \
        --pred regenie_step1_chr${chrom}_pred.list \
        --chr ${chrom} \
        --out regenie_chr${chrom}.out
    """
}

// 🔹 Workflow principal
workflow {
    qc_results = quality_control(channel_files)
    regenie_results = run_regenie(qc_results.bed, qc_results.bim, qc_results.fam, channel_chroms)

    regenie_results.view() // 📢 Affiche les fichiers générés
}
