meta="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/data/CCGP_Metadata_Submission_C_fasciata_Nachman-Bowie_v2_oct2022_MAT.csv"
outdir="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/ccgp_chamaea/scripts"

mkdir -p "$outdir"

awk -F',' -v outdir="$outdir" '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "sampleID") sample_col = i
        if ($i == "Ecoregion") eco_col = i
    }
    next
}

{
    sample = $sample_col
    eco = $eco_col

    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", sample)
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", eco)

    if (sample == "") next
    if (eco == "") eco = "Unknown"

    n = split(eco, words, /[^A-Za-z0-9]+/)
    acronym = ""

    for (i = 1; i <= n; i++) {
        if (words[i] != "") {
            acronym = acronym toupper(substr(words[i], 1, 1))
        }
    }

    if (acronym == "") acronym = "UNKNOWN"

    file = outdir "/" acronym ".txt"
    print sample >> file
}
' "$meta"