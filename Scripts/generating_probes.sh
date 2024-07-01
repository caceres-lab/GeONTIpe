#!/bin/bash

wd="$(pwd)"

region=$1
extra_region=$2
reference=$3
reference2=$4
liftOver_route=$5

mkdir -p $wd/Infor/blast_results

LIFTOVER_CHAIN="$liftOver_route/hg38tot2t.over.chain"

blast_window() {
    local window_file=$1
    
    echo $window_file
    
    region=$(echo $window_file | cut -d'_' -f2)
    type=$(echo $window_file | cut -d'_' -f3)
    chr=$(echo $window_file | cut -d'_' -f4)

    region_start="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f3)"
    region_end="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f6)"
    
    result_file="$wd/Infor/blast_results/$window_file.out"

    endchr="$(cat $wd/Infor/conversion.txt | grep $chr[[:space:]] | cut -f3)"
    newchrname="$(cat $wd/Infor/conversion.txt | grep $chr[[:space:]] | cut -f1)"

    region_start_ext=$((region_start - $extra_region))
    region_end_ext=$((region_end + $extra_region))
    
    if [ $region_start_ext -lt 1 ]; then
        region_start_ext=1
    fi

    if [ $region_end_ext -gt $endchr ]; then
        region_end_ext=$endchr
    fi

    seq_region=$(samtools faidx $reference "$newchrname:${region_start_ext}-${region_end_ext}")
    echo "$seq_region" > $wd/Infor/windows/$region/region_tmp.fasta

    if [ -s $wd/Infor/windows/$region/$type/$window_file ] || [ -s $wd/Infor/windows/$region/region_tmp.fasta ];then
      blastn -task megablast -qcov_hsp_perc 70 -query $wd/Infor/windows/$region/$type/$window_file -subject $wd/Infor/windows/$region/region_tmp.fasta -outfmt 6 > $result_file
    else 
      > $result_file
    fi
    
    hits=$(cat $result_file | wc -l)
    
    if (( hits == 1 )); then
        echo "$window_file hg38_correct_probe" >> $wd/Infor/windows/$region/tracking_windows.txt
        
        result_file_t2t="$wd/Infor/blast_results/${window_file}_t2t.out"

        echo -e "${chr}\t$region_start\t$region_end" > $wd/Infor/windows/$region/coords.bed
        name_region_conversion="$(echo "${chr}_${region_start}_${region_end}")"
        
        if [ ! -f $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed ];then
            liftOver -minMatch=0.25 -multiple $wd/Infor/windows/$region/coords.bed $LIFTOVER_CHAIN $wd/Infor/windows/$region/coords_t2t.bed $wd/Infor/windows/$region/Unmapped_t2t.bed
        fi
            
        if [ -s $wd/Infor/windows/$region/coords_t2t.bed ];then
          chr="$(head -1 $wd/Infor/windows/$region/coords_t2t.bed | cut -f1)"
          st="$(head -1 $wd/Infor/windows/$region/coords_t2t.bed | cut -f2)"
          newst=$((st - extra_region))
          end="$(head -1 $wd/Infor/windows/$region/coords_t2t.bed | cut -f3)"
          newend=$((end + extra_region))
          echo -e "$chr\t$newst\t$newend" > $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed
        
          t2t_region=$(awk '{print $1 ":" $2 "-" $3}' $wd/Infor/windows/$region/${name_region_conversion}_coords_t2t.bed)
          result_file_t2t="$wd/Infor/blast_results/${window_file}_t2t.out"
          seq_region_t2t=$(samtools faidx $reference2 "$t2t_region")
          echo "$seq_region_t2t" > $wd/Infor/windows/$region/region_tmp_t2t.fasta
          
          if [ -s $wd/Infor/windows/$region/$type/$window_file ] || [ -s $wd/Infor/windows/$region/region_tmp_t2t.fasta ];then
            blastn -task megablast -qcov_hsp_perc 70 -query $wd/Infor/windows/$region/$type/$window_file -subject $wd/Infor/windows/$region/region_tmp_t2t.fasta -outfmt 6 > $result_file_t2t
          else 
            > $result_file
          fi
          
          hits_t2t=$(cat $result_file_t2t | wc -l)
          
          if (( hits_t2t == 1 )); then
            echo $window_file > $wd/Infor/windows/$region/$type/specific_windows.txt
            echo "$window_file t2t_correct_probe" >> $wd/Infor/windows/$region/tracking_windows.txt
            return 0 
          else
            echo "$window_file t2t_wrong_probe" >> $wd/Infor/windows/$region/tracking_windows.txt
            rm $wd/Infor/windows/$region/$type/$window_file
          fi
        else
          echo "$window_file wrong_hg38_t2t_conversion" >> $wd/Infor/windows/$region/tracking_windows.txt
          rm $wd/Infor/windows/$region/$type/$window_file
        fi
    else
      echo "$window_file hg38_wrong_probe" >> $wd/Infor/windows/$region/tracking_windows.txt
      rm $wd/Infor/windows/$region/$type/$window_file
    fi

    return 1
}

export -f blast_window

found_probe=false

if [ ! -s $wd/Infor/windows/$region/$type/allwindows.txt ];then
  echo "$region $type not windows" >> $wd/Infor/probecreation/generated_probes.txt
else
  while read -r probe; do
      region_info=$(basename $probe)
        if blast_window "$region_info"; then
          found_probe=true
          break
        fi
  done < $wd/Infor/windows/$region/allwindows.txt
fi

if $found_probe; then
    find $wd/Infor/windows/$region/$type -type f ! -name "$(head -n 1 $wd/Infor/windows/$region/$type/specific_windows.txt | xargs basename)" -delete
    
    mkdir -p $wd/Genotyping
    mkdir -p $wd/Genotyping/$region
    mkdir -p $wd/Genotyping/$region/CoordenadasRelativas
    mkdir -p $wd/Genotyping/$region/Sondas
    
    if [ $(ls $wd/Infor/windows/$region/$type | wc -l) == 0 ];then
        echo "$region $type No probes" >> $wd/Infor/probecreation/generated_probes.txt
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
        echo "$region $type found probes" >> $wd/Infor/probecreation/generated_probes.txt
    fi
else
  echo "$region $type not windows" >> $wd/Infor/probecreation/generated_probes.txt
fi
