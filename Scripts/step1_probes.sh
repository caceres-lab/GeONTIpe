#!/bin/bash

wd="$(pwd)"

flanking_region=$1
masked_genome_route=$2
liftOver_route=$2

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
