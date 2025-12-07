# Make multiple "gene" trees per chromosome

ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas/
NC_087513.1.fa  NC_087517.1.fa  NC_087521.1.fa  NC_087525.1.fa  NC_087529.1.fa  NC_087533.1.fa  NC_087537.1.fa  NC_087541.1.fa  NC_087545.1.fa
NC_087514.1.fa  NC_087518.1.fa  NC_087522.1.fa  NC_087526.1.fa  NC_087530.1.fa  NC_087534.1.fa  NC_087538.1.fa  NC_087542.1.fa  NC_087546.1.fa
NC_087515.1.fa  NC_087519.1.fa  NC_087523.1.fa  NC_087527.1.fa  NC_087531.1.fa  NC_087535.1.fa  NC_087539.1.fa  NC_087543.1.fa
NC_087516.1.fa  NC_087520.1.fa  NC_087524.1.fa  NC_087528.1.fa  NC_087532.1.fa  NC_087536.1.fa  NC_087540.1.fa  NC_087544.1.fa


#!/bin/bash

# Set your input and output directories
input_dir="/xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas"
output_dir="/xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas/raxml_windows"

mkdir -p "$output_dir"

for fasta in "$input_dir"/*.fa; do
    chrom=$(basename "$fasta" .fa)
    echo "Processing $chrom"

    # Get sequence length
    seqlen=$(awk '/^[^>]/ {seqlen += length($0)} END {print seqlen}' "$fasta")

    # Split into 50kb windows
    start=1
    window_size=50000

    while [ $start -le $seqlen ]; do
        end=$((start + window_size - 1))
        if [ $end -gt $seqlen ]; then
            end=$seqlen
        fi

        # Extract window
        output="$output_dir/${chrom}_${start}_${end}.fa"
        seqkit subseq -r $start:$end "$fasta" > "$output"

        echo "Created $output"

        start=$((end + 1))
    done
done
