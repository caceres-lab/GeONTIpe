#!/bin/bash

wd="$(pwd)"

REGIONS_FILE=$1
WINDOW_SIZE=$2
OVERLAP=$3
STEP=$((WINDOW_SIZE - OVERLAP))
masked_genome_route=$4
minim_rep_seq=$5

masked_genome="$masked_genome_route/hg38_masked.fa"

create_windows() {
    local seq=$1
    local region_name=$2
    local region_start=$3
    local length=${#seq}
    for ((i=0; i<length-WINDOW_SIZE+1; i+=STEP)); do
        local window_seq=${seq:i:WINDOW_SIZE}
        local start=$((region_start + i))
        local end=$((start + WINDOW_SIZE - 1))
        local window_name="${region_name}_${start}_${end}"
        local inversion=$(echo $region_name | cut -d'_' -f1)
        local type=$(echo $region_name | cut -d'_' -f2)
        mkdir -p $wd/Infor/windows/$inversion/$type

        local lower_count=$(echo $window_seq | grep -o '[acgt]' | wc -l)
        local total_count=${#window_seq}
        local lower_percentage=$((100 * lower_count / total_count))

        if [[ $lower_percentage -lt $minim_rep_seq ]] && [[ $window_seq != *"N"* ]]; then
            echo ">${window_name}" > $wd/Infor/windows/$inversion/$type/window_${window_name}.fasta
            echo $window_seq >> $wd/Infor/windows/$inversion/$type/window_${window_name}.fasta
        fi
    done
    
    echo "$region ready" >> $wd/tracking_pipeline/windows_ready.txt
}

region=$(echo $REGIONS_FILE | awk '{print $1}')
chr=$(echo $REGIONS_FILE | awk '{print $2}')
region_coords=$(echo $REGIONS_FILE | awk '{print $3}')
region_start=$(echo $region_coords | cut -d'-' -f1)
seq=$(samtools faidx $masked_genome "chr${chr}:${region_coords}" | sed -n '1!p' | tr -d '\n')

create_windows "$seq" "${region}_chr${chr}" $region_start
