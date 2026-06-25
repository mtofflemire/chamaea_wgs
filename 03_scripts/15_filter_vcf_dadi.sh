#filter out the first scaffold (chrom1)
bcftools view \
  -r JARCOQ010000001.1 \
  -O z \
  -o /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_raw_chr1.vcf.gz \
  /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/ccgp/qc/41-Chamaea_raw.vcf.gz

#just do quality control filtering
bcftools view \
  -v snps \
  -f .,PASS \
  -e 'ALT="*" | TYPE~"indel" | ref="N" | QUAL < 30.0 | MQ < 40.0 | MQRankSum < -12.5 | FS > 60.0 | SOR > 3.0 | ReadPosRankSum < -8.0 | QD < 2.0' \
  -O z \
  -o /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_filteredQC_gatk4_chr1.vcf.gz \
  /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_raw_chr1.vcf.gz

#filters for dadi analyis
bcftools view \
  -S /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea-samples.txt \
  -t ^mtDNA \
  -v snps \
  -m2 -M2 \
  -f .,PASS \
  -e 'AF==1 | AF==0 | ALT="*" | F_MISSING > 0 | TYPE~"indel" | ref="N" | QUAL < 30.0 | MQ < 40.0 | MQRankSum < -12.5 | FS > 60.0 | SOR > 3.0 | ReadPosRankSum < -8.0 | QD < 2.0' \
  -O z \
  -o /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_filteredQC_chr1.vcf.gz \
  /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_raw_chr1.vcf.gz


#check to see how many snps are in your QC_SNPs
bcftools view -H /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_filteredQC_gatk4_chr1.vcf.gz | wc -l

#COUNT SNPS: check to see how many snps are in your filtered_QC_SNPs
bcftools view -H /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_raw_chr1.vcf.gz | wc -l


#Calculate length (L) of chrom1=total_callable_sites*(filtered_QC_SNPs/total_QC_SNPS)
#L=151832209*(4804850/5182481)
#L=140768676.125
#Ne= 883331.4848443555/(4*(2.3e-9)*140,768,676.13)
#Ne= 682071.427185


#best model for Secondary Contact
#Ne= 883331.4848443555/(4*(2.3e-9)*140,768,676.13)
#Ne=682071.427185
#Ne-north=4.29775664e-01*682071.427185 = 248808.978923
#Ne-south=0.99269426868*682071.427185 = 750946.810711
#T-no-mig =0.17880423*2*307720.594247*1 = 110043.487819
#T-mig=4.02633479e-01*2*682071.427185*1 = 82597.7218489



Step 1: Compute Migration Rates
For  m12 :
m12 = \frac{3.5250106}{2 \times 307720.594247}
m12 = \frac{3.5250106}{615441.188494} = 5.727 \times 10^{-6}
For  m21 :
m21 = \frac{3.89418467}{2 \times 307720.594247}
m21 = \frac{3.89418467}{615441.188494} = 6.327 \times 10^{-6}
tep 2: Compute Absolute Number of Migrants

For south → north (Pop2 → Pop1):
N_{m12} = m12 \times N_{\text{dest}} = (5.727 \times 10^{-6}) \times 248808.978923
N_{m12} = 1.424
For north → south (Pop1 → Pop2):
N_{m21} = m21 \times N_{\text{dest}} = (6.327 \times 10^{-6}) \times 750946.810711
N_{m21} = 4.753
