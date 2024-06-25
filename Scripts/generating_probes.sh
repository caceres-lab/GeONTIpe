#!/bin/bash

wd="$(pwd)"

flanking_region=$1
extra_region=$2
WINDOW_SIZE=$3
OVERLAP=$4
STEP=$((WINDOW_SIZE - OVERLAP))
reference=$5
reference2=$6
masked_genome_route=$7
liftOver_route=$7
minim_rep_seq=$8

mkdir $wd/Infor/probecreation

if [ ! -f $wd/Infor/Reference/hg38_masked.fa ];then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip -c hg38.fa.gz > $masked_genome_route/hg38_masked.fa
    rm hg38.fa.gz
fi

if [ ! -f $wd/Infor/Reference/hg38tot2t.over.chain ];then
    wget https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz
    gunzip -c hg38-chm13v2.over.chain.gz > $liftOver_route/hg38tot2t.over.chain
    rm hg38-chm13v2.over.chain.gz
fi

masked_genome="$masked_genome_route/hg38_masked.fa"
LIFTOVER_CHAIN="$liftOver_route/hg38tot2t.over.chain"

merging_BC=$((flanking_region*2))

> $wd/Infor/probecreation/regions_probes.txt

for region in $(cut -f1 $wd/Infor/allcoords.txt)
do
    if [ -d $wd/Genotyping/$region ] && [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 4 ];then
        continue
    else
        chr="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f2)"
        invsize="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f7)"

        if [ $chr == 23 ];then
            chr="X"
        elif [ $chr == 24 ];then
            chr="Y"
        fi

        endchr="$(cat $wd/Infor/conversion.txt | grep chr$chr[[:space:]] | cut -f3)"
        endA="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f3)"
        startA="$((endA - flanking_region))"

        if [ $startA -le 0 ];then
            startA=1
        fi

        regionA="$startA-$endA"
        echo -e "${region}_A\t$chr\t$regionA" >> $wd/Infor/probecreation/regions_probes.txt

        if [ $invsize -le $merging_BC ];then
            startB="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f4)"
            newinv_size=$((invsize/2))
            endB="$((startB + newinv_size))"
            regionB="$startB-$endB"
            
            endC="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f5)"
            startC="$((endC - newinv_size))"
            regionC="$startC-$endC"
            
            echo -e "${region}_B\t$chr\t$regionB" >> $wd/Infor/probecreation/regions_probes.txt
            echo -e "${region}_C\t$chr\t$regionC" >> $wd/Infor/probecreation/regions_probes.txt
        else
            startB="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f4)"
            endB="$((startB + flanking_region))"

            if [ $endB -ge $endchr ];then
                endB=$endchr
            fi

            regionB="$startB-$endB"

            endC="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f5)"
            startC="$((endC - flanking_region))"

            regionC="$startC-$endC"
            echo -e "${region}_B\t$chr\t$regionB" >> $wd/Infor/probecreation/regions_probes.txt
            echo -e "${region}_C\t$chr\t$regionC" >> $wd/Infor/probecreation/regions_probes.txt
        fi
        
        startD="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f6)"
        endD="$((startD + flanking_region))"

        if [ $endD -ge $endchr ];then
            endD=$endchr
        fi

        regionD="$startD-$endD"
        echo -e "${region}_D\t$chr\t$regionD" >> $wd/Infor/probecreation/regions_probes.txt
        echo -e "${region}\t$startA" >> $wd/Infor/probecreation/beginning_region.txt
    fi
done

REGIONS_FILE="$wd/Infor/probecreation/regions_probes.txt"

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
}

while read -r line; do
    region=$(echo $line | awk '{print $1}')
    chr=$(echo $line | awk '{print $2}')
    region_coords=$(echo $line | awk '{print $3}')
    region_start=$(echo $region_coords | cut -d'-' -f1)
    echo $region
    seq=$(samtools faidx $masked_genome "chr${chr}:${region_coords}" | sed -n '1!p' | tr -d '\n')

    create_windows "$seq" "${region}_chr${chr}" $region_start
done < $REGIONS_FILE

for region in $(cut -f1 $wd/Infor/allcoords.txt)
do
    if [ -d $wd/Genotyping/$region ] && [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 4 ];then
        continue
    else
      if [ -f $wd/Infor/windows/$region/A ];then
        ls $wd/Infor/windows/$region/A | grep -v "allwindows" | sort -rn > $wd/Infor/windows/$region/A/allwindows.txt
        echo $region/A >> $wd/Infor/probecreation/allregions.txt
      fi
      
      if [ -f $wd/Infor/windows/$region/B ];then
        ls $wd/Infor/windows/$region/B | grep -v "allwindows" | sort -n > $wd/Infor/windows/$region/B/allwindows.txt
        echo $region/B >> $wd/Infor/probecreation/allregions.txt
      fi
      
      if [ -f $wd/Infor/windows/$region/C ];then
        ls $wd/Infor/windows/$region/C | grep -v "allwindows" | sort -rn > $wd/Infor/windows/$region/C/allwindows.txt
        echo $region/C >> $wd/Infor/probecreation/allregions.txt
      fi
      
      if [ -f $wd/Infor/windows/$region/D ];then
        ls $wd/Infor/windows/$region/D | grep -v "allwindows" | sort -n > $wd/Infor/windows/$region/D/allwindows.txt
        echo $region/D >> $wd/Infor/probecreation/allregions.txt
      fi
    fi
done

makeblastdb -in $reference -dbtype nucl
makeblastdb -in $reference2 -dbtype nucl

mkdir -p $wd/Infor/blast_results

LIFTOVER_CHAIN="$wd/Infor/Reference/hg38tot2t.over.chain"

blast_window() {
    local window_file=$1
    
    region=$(echo $window_file | cut -d'_' -f2)
    type=$(echo $window_file | cut -d'_' -f3)
    chr=$(echo $window_file | cut -d'_' -f4)

    region_start="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f3)"
    region_end="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f6)"
    
    result_file="$wd/Infor/blast_results/$window_file.out"

    endchr="$(cat $wd/Infor/conversion.txt | grep $chr[[:space:]] | cut -f3)"
    newchrname="$(cat $wd/Infor/conversion.txt | grep $chr[[:space:]] | cut -f2)"

    region_start_ext=$((region_start - $extra_region))
    region_end_ext=$((region_end + $extra_region))
    
    if [ $region_start_ext -lt 1 ]; then
        region_start_ext=1
    fi

    if [ $region_end_ext -gt $endchr ]; then
        region_end_ext=$endchr
    fi

    seq_region=$(samtools faidx $reference "$newchrname:${region_start_ext}-${region_end_ext}")
    echo "$seq_region" > $wd/region_tmp.fasta

    blastn -query $wd/Infor/windows/$region/$type/$window_file -subject $wd/region_tmp.fasta -outfmt 6 > $result_file

    hits=$(cat $result_file | wc -l)
    if (( hits == 1 )); then
        result_file_t2t="$wd/Infor/blast_results/${window_file}_t2t.out"

        echo -e "${chr}\t$region_start\t$region_end" > $wd/Infor/coords.bed
        name_region_conversion="$(echo "${chr}_${region_start}_${region_end}")"
        
        if [ ! -f $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed ];then
            liftOver -minMatch=0.8 -multiple $wd/Infor/coords.bed $LIFTOVER_CHAIN $wd/Infor/windows/$region/coords_t2t.bed $wd/Infor/windows/$region/Unmapped_t2t.bed
            chr="$(head -1 $wd/Infor/windows/$region/coords_t2t.bed | cut -f1)"
            st="$(head -1 $wd/Infor/windows/$region/coords_t2t.bed | cut -f2)"
            newst=$((st - extra_region))
            end="$(head -1 $wd/Infor/windows/$region/coords_t2t.bed | cut -f3)"
            newend=$((end + extra_region))
            echo -e "$chr\t$newst\t$newend" > $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed
        fi

        if [ -s $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed ] && [ -s $wd/Infor/windows/$region/coords_t2t.bed ]; then
            t2t_region=$(awk '{print $1 ":" $2 "-" $3}' $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed)
            result_file_t2t="$wd/Infor/blast_results/${window_file}_t2t.out"
            seq_region_t2t=$(samtools faidx $reference2 "$t2t_region")
            echo "$seq_region_t2t" > $wd/region_tmp_t2t.fasta
            
            blastn -query $wd/Infor/windows/$region/$type/$window_file -subject $wd/region_tmp_t2t.fasta -outfmt 6 > $result_file_t2t

            hits_t2t=$(cat $result_file_t2t | wc -l)
            if (( hits_t2t == 1 )); then
                echo $window_file > $wd/Infor/windows/$region/$type/specific_windows.txt
                return 0 
            fi
        fi
    fi

    return 1
}

export -f blast_window

while read -r reg; do
    found_probe=false
    while read -r probe; do
        region_info=$(basename $probe)
        if blast_window "$region_info"; then
            found_probe=true
            break
        fi
    done < $wd/Infor/windows/$reg/allwindows.txt

    if $found_probe; then
        find $wd/Infor/windows/$reg -type f ! -name "$(head -n 1 $wd/Infor/windows/$region/$type/specific_windows.txt | xargs basename)" -delete
        mkdir -p $wd/Genotyping
        mkdir -p $wd/Genotyping/$region
        mkdir -p $wd/Genotyping/$region/CoordenadasRelativas
        mkdir -p $wd/Genotyping/$region/Sondas
        if [ $(ls $wd/Infor/windows/$region/$type | wc -l) == 0 ];then
            echo "$region $type No probe" >> $wd/Infor/probecreation/generated_probes.txt
        elif [ $(ls $wd/Infor/windows/$region/$type | wc -l) -ge 2 ];then
            echo "$region $type Much probes" >> $wd/Infor/probecreation/generated_probes.txt
        else
            cat $wd/Infor/windows/$region/$type/* | head -2 | tail -1 > $wd/Genotyping/$region/Sondas/$type
            beg="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
            coord1="$(echo $wd/Infor/windows/$region/$type/* | sed 's/_/\t/g' | cut -f5)"
            coord2="$(echo $wd/Infor/windows/$region/$type/* | sed 's/_/\t/g' | cut -f6 | sed 's/.fasta/ /g')"
            relativecoord1=$((coord1 - beg))
            relativecoord2=$((coord2 - beg))
            echo -e "$relativecoord1\t$relativecoord2" > $wd/Genotyping/$region/CoordenadasRelativas/$type
        fi
    fi
done < $wd/Infor/probecreation/allregions.txt

rm *.fasta

> $wd/tracking_pipeline/created_probes.txt

for region in $(cut -f1 $wd/Infor/allcoords.txt)
do
    if [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 4 ];then
        echo "$region Probes completed" >> $wd/tracking_pipeline/created_probes.txt
    elif [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 3 ];then
        echo "$region Set partial" >> $wd/tracking_pipeline/created_probes.txt
        if [ ! -f $wd/Genotyping/$region/A ];then
            bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
            co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f3)"
            relcoorda1=$((co1 - bega - 500))
            relcoorda2=$((relcoorda1 + 300))
            echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/A
            > $wd/Genotyping/$region/Sondas/A
        elif [ ! -f $wd/Genotyping/$region/B ];then
            bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
            co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f4)"
            relcoorda1=$((co1 - bega + 500))
            relcoorda2=$((relcoorda1 + 300))
            echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/B
            > $wd/Genotyping/$region/Sondas/B
        elif [ ! -f $wd/Genotyping/$region/C ];then
            bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
            co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f5)"
            relcoorda1=$((co1 - bega - 500))
            relcoorda2=$((relcoorda1 + 300))
            echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/C
            > $wd/Genotyping/$region/Sondas/C
        elif [ ! -f $wd/Genotyping/$region/D ];then
            bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
            co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f6)"
            relcoorda1=$((co1 - bega + 500))
            relcoorda2=$((relcoorda1 + 300))
            echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/D
            > $wd/Genotyping/$region/Sondas/D
        fi
    else 
        echo "$region Failed probes" >> $wd/tracking_pipeline/created_probes.txt
        rm -r $wd/Genotyping/$region
    fi
done

#rm -r $wd/Infor/probecreation
#rm -r $wd/Infor/windows
#rm -r $wd/Infor/blast_results