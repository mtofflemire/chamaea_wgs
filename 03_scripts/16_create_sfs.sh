#!/bin/bash

# Define environment name and desired packages
ENV_NAME="easySFS_env"
PYTHON_VERSION="3.9"
PACKAGES="scipy numpy pandas"

# Initialize Conda (adjust the path if your conda is located elsewhere)
source ~/anaconda3/etc/profile.d/conda.sh

# Check if the environment already exists
if conda env list | grep -qE "^$ENV_NAME\s"; then
    echo "Environment '$ENV_NAME' already exists. Activating it..."
else
    echo "Environment '$ENV_NAME' does not exist. Creating it now..."
    conda create -n "$ENV_NAME" python=$PYTHON_VERSION $PACKAGES -y
fi

# Activate the environment
conda activate "$ENV_NAME"

# Verify dependencies
python -c "import scipy; import numpy; import pandas; print('All modules loaded successfully!')"

# Define paths
VCF="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_filteredQC_chr1.vcf.gz"
POPFILE="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea.pop.clade_north_south_copy.txt"
OUTDIR="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results"

# Run easySFS in preview mode
echo "Running easySFS in preview mode to suggest projection sizes..."
python /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/systems/easySFS/easySFS.py \
  -i "$VCF" \
  -p "$POPFILE" \
  -o "$OUTDIR" \
  --preview -a  # Include all SNPs

# Prompt the user for projection sizes
echo "------------------------------------------------------------"
echo "Preview complete. Please enter the projection sizes."
echo "Format: <projection_pop1>,<projection_pop2> (e.g., 40,40)"
read -p "Enter projection sizes: " PROJ_SIZES

# Confirm projection size selection
echo "------------------------------------------------------------"
echo "Using projection sizes: $PROJ_SIZES"
read -p "Press Enter to continue, or Ctrl+C to cancel..."

# Run easySFS with selected projections
python /Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/systems/easySFS/easySFS.py \
  -i "$VCF" \
  -p "$POPFILE" \
  -o "$OUTDIR" \
  --proj "$PROJ_SIZES" -a  # Keep all SNPs

echo "------------------------------------------------------------"
echo "SFS generation complete! Files saved in: $OUTDIR"





#this is a quick way of inpurting your particular project size based on code from easySFS website
#VCF="/Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_filteredQC_chr1.vcf.gz"
#POPFILE="/Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea.pop.clade_north_south_copy.txt"
#OUTDIR="/Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/results"#

#preview projection sizes for maximizing number of segregating sites for each population
#python /Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/systems/easySFS/easySFS.py \
#  -i "$VCF" \
#  -p "$POPFILE" \
#  -o "$OUTDIR" \
#  --preview -a  


#this is a quick way of inpurting your particular project size based on code from easySFS website
#VCF="/Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea_filteredQC_chr1.vcf.gz"
#POPFILE="/Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea.pop.clade_north_south_copy.txt"
#OUTDIR="/Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/results"#

#python /Users/mtofflemire/Dropbox/sites/storage/local/projects/chamaea_genomics/systems/easySFS/easySFS.py \
#  -i "$VCF" \
#  -p "$POPFILE" \
#  -o "$OUTDIR" \
#  --proj 168,150 -a  