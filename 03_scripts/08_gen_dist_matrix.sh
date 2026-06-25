#This is the one that worked for the EEMS analysis. For some reason, problems arose from the ones above. 
plink --bfile /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/Chamaea_genomics/data/41-Chamaea_filteredQC_autosomal_only \
      --distance square 1-ibs \
      --out Chamaea_filteredQC_autosomal_eems
