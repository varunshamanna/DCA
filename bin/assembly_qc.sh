# Extract assembly QC information and determine QC result based on report.tsv from Quast, total base count
BASES=$(< $JSON jq -r .summary.after_filtering.total_bases)
LENGTH=$(awk -F'\t' '$1 == "Total length" { print $2 }' $REPORT)
ASSEMBLY_DEPTH=$(printf %.2f $(echo "$BASES / $LENGTH" | bc -l) )

if (( $(echo "$ASSEMBLY_DEPTH >= $QC_DEPTH" | bc -l) )); then
    ASSEMBLY_DEPTH_QC="PASS"
else
    ASSEMBLY_DEPTH_QC="FAIL"
fi
