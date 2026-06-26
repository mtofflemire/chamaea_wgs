#!/bin/bash
set -euo pipefail



#define pathways needed for pca
SNP_DIR='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/2_data/snps'
OUTPCA='/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs/4_analyses/01_pca'
BCFTOOLS='/Users/michaeltofflemire/miniconda3/envs/vcf/bin/bcftools'
PLINK='/Users/michaeltofflemire/miniconda3/envs/vcf/bin/plink'
mkdir -p "$SNP_DIR"
mkdir -p "$OUTPCA"


#define qc filtered vcf file and prefix
INVCF=$SNP_DIR/Chamaea_auto_filteredQC.vcf.gz
VCF=$SNP_DIR/Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb.vcf.gz
BED=$OUTPCA/Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb
PCA=$OUTPCA/Chamaea_auto_filteredQC_maf0.05_hwe0.01_prune1kb.pca



#add population filters for autosomal structure analysis
"$BCFTOOLS" +fill-tags "$INVCF" -- -t MAF,HWE \
| "$BCFTOOLS" view -t ^CM055129.1,JARCOQ010000004.1,JARCOQ010000029.1 \
| "$BCFTOOLS" view -v snps -m2 -M2 \
  -i 'MAF[0]>=0.05 && HWE>=0.01' \
| "$BCFTOOLS" +prune -n 1 -N rand --random-seed 12345 -w 1kb -Oz -o "$VCF"



#index vcf
"$BCFTOOLS" index --force "$VCF"



#run pca with PLINK
"$PLINK" --vcf "$VCF" --double-id --allow-extra-chr --make-bed --threads 8 --out "$BED"
"$PLINK" --bfile "$BED" --allow-extra-chr --pca 20 header tabs --threads 8 --out "$PCA"
