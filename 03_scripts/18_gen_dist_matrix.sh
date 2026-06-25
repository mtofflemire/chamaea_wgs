cd "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea"

plink \
  --vcf "results/01_VARIANTS/GATK_QC/AUTOSOMES/Chamaea_autoALL_filteredQC_maf0.05_miss0.25_prune1kb.vcf.gz" \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out "scripts/chamaea_plink"

plink \
  --bfile "scripts/chamaea_plink" \
  --allow-extra-chr \
  --distance square 1-ibs \
  --out "scripts/chamaea_gendist"