# Return reference coverage percentage by the reads

samtools index "$BAM"
COVERAGE=$(samtools coverage "$BAM" | awk -F'\t' 'FNR==2 {print $6}')
DEPTH=$(samtools depth -a "$BAM" | awk '{c++;s+=$3}END{print s/c}')