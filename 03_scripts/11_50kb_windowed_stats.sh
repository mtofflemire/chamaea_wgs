#!/bin/bash

# Define paths
VCF_FILE="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea-filteredQC.vcf.gz"
POPULATION_FILE="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea.pop.clade_north_south_copy.txt"
DATA_DIR="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data"
RESULTS_DIR="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/10kb-windows-results"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Define paths to permanent population files (stored in data/)
POP_SOUTH_FILE="$DATA_DIR/41-Chamaea.pops_South.txt"
POP_NORTH_FILE="$DATA_DIR/41-Chamaea.pops_North.txt"

# Define output prefixes
FST_OUTPUT_PREFIX="$RESULTS_DIR/south_v_north_fst"
PI_OUTPUT_PREFIX="$RESULTS_DIR/south_v_north_pi"
TAJIMA_OUTPUT_PREFIX="$RESULTS_DIR/south_v_north_tajimasD"

### Step 1: Check for Existing Population Files Before Extracting ###
if [[ -s "$POP_SOUTH_FILE" && -s "$POP_NORTH_FILE" ]]; then
    echo "✅ Population files already exist in $DATA_DIR. Skipping extraction."
else
    echo "🔹 Extracting permanent population files..."
    
    # Extract individuals into permanent population files
    awk '$2 == "South" {print $1}' "$POPULATION_FILE" > "$POP_SOUTH_FILE"
    awk '$2 == "North" {print $1}' "$POPULATION_FILE" > "$POP_NORTH_FILE"

    # Debugging: Check if files were created successfully
    if [[ ! -s "$POP_SOUTH_FILE" || ! -s "$POP_NORTH_FILE" ]]; then
        echo "❌ ERROR: One or both population files are empty or missing!"
        echo "South population count: $(wc -l < "$POP_SOUTH_FILE" 2>/dev/null || echo 0)"
        echo "North population count: $(wc -l < "$POP_NORTH_FILE" 2>/dev/null || echo 0)"
        exit 1
    fi

    echo "✅ South population count: $(wc -l < "$POP_SOUTH_FILE")"
    echo "✅ North population count: $(wc -l < "$POP_NORTH_FILE")"
fi

### Step 2: Run vcftools Fst calculation ###
#echo "🔹 Calculating Fst..."
#vcftools --gzvcf "$VCF_FILE" \
#         --weir-fst-pop "$POP_SOUTH_FILE" \
#         --weir-fst-pop "$POP_NORTH_FILE" \
#         --fst-window-size 10000 \
#         --fst-window-step 10000 \
#         --out "$FST_OUTPUT_PREFIX"

#echo "✅ Fst calculation complete. Output saved in $FST_OUTPUT_PREFIX"

### Step 3: Run vcftools windowed π calculation ###
echo "🔹 Calculating windowed π (nucleotide diversity)..."

# Calculate windowed π for the South population
#vcftools --gzvcf "$VCF_FILE" \
#         --keep "$POP_SOUTH_FILE" \
#         --window-pi 10000 \
#         --out "${PI_OUTPUT_PREFIX}_south"

# Calculate windowed π for the North population
#vcftools --gzvcf "$VCF_FILE" \
#         --keep "$POP_NORTH_FILE" \
#         --window-pi 10000 \
#         --out "${PI_OUTPUT_PREFIX}_north"

#echo "✅ Windowed π calculation complete. Outputs saved in:"
#echo "  - ${PI_OUTPUT_PREFIX}_south"
#echo "  - ${PI_OUTPUT_PREFIX}_north"

### Step 4: Run vcftools Tajima’s D calculation ###
echo "🔹 Calculating Tajima’s D..."

# Calculate Tajima’s D for the South population
vcftools --gzvcf "$VCF_FILE" \
         --keep "$POP_SOUTH_FILE" \
         --TajimaD 10000 \
         --out "${TAJIMA_OUTPUT_PREFIX}_south"

# Calculate Tajima’s D for the North population
vcftools --gzvcf "$VCF_FILE" \
         --keep "$POP_NORTH_FILE" \
         --TajimaD 10000 \
         --out "${TAJIMA_OUTPUT_PREFIX}_north"

echo "✅ Tajima’s D calculation complete. Outputs saved in:"
echo "  - ${TAJIMA_OUTPUT_PREFIX}_south.Tajima.D"
echo "  - ${TAJIMA_OUTPUT_PREFIX}_north.Tajima.D"

echo "🎉 All analyses finished successfully!"







#Use this to copy and paste into bash for simple tajima's D calculate for a particular population
#!/bin/bash

# Define input files
VCF_FILE="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/41-Chamaea-filteredQC.vcf.gz"       # Replace with your VCF file path
POP_FILE="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/data/fst_sample_lists/oregon_samples.txt"    # Replace with your population file path
OUTPUT_PREFIX="/Users/michaeltofflemire/Library/CloudStorage/Dropbox/sites/storage/local/projects/chamaea_genomics/results/chamaea-10kb-windowed/Oregon_tajimasD"  # Replace with desired output prefix

# Ensure output directory exists
mkdir -p "$(dirname "$OUTPUT_PREFIX")"

# Run Tajima's D calculation
vcftools --gzvcf "$VCF_FILE" \
         --keep "$POP_FILE" \
         --TajimaD 10000 \
         --out "$OUTPUT_PREFIX"

echo "✅ Tajima’s D calculation complete. Output saved in:"
echo "  - ${OUTPUT_PREFIX}.Tajima.D"
