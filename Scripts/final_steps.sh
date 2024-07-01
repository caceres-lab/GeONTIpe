#!/bin/bash

wd="$(pwd)"

> $wd/tracking_pipeline/created_probes.txt

for region in $(cut -f1 $wd/Infor/allcoords.txt)
do
    
    if [ ! -d $wd/Genotyping/$region ];then
      echo "$region Failed probes" >> $wd/tracking_pipeline/created_probes.txt
    else
      if [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 4 ];then
          echo "$region Probes completed" >> $wd/tracking_pipeline/created_probes.txt
      elif [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 3 ];then
          echo "$region Set partial" >> $wd/tracking_pipeline/created_probes.txt
          if [ ! -f $wd/Genotyping/$region/Sondas/A ];then
              bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
              co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f3)"
              relcoorda1=$((co1 - bega - 500))
              relcoorda2=$((relcoorda1 + 300))
              echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/A
              > $wd/Genotyping/$region/Sondas/A
          elif [ ! -f $wd/Genotyping/$region/Sondas/B ];then
              bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
              co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f4)"
              relcoorda1=$((co1 - bega + 500))
              relcoorda2=$((relcoorda1 + 300))
              echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/B
              > $wd/Genotyping/$region/Sondas/B
          elif [ ! -f $wd/Genotyping/$region/Sondas/C ];then
              bega="$(cat $wd/Infor/probecreation/beginning_region.txt | grep $region[[:space:]] | cut -f2)"
              co1="$(cat $wd/Infor/allcoords.txt | grep ${region}[[:space:]] | cut -f5)"
              relcoorda1=$((co1 - bega - 500))
              relcoorda2=$((relcoorda1 + 300))
              echo -e "$relcoorda1\t$relcoorda2" > $wd/Genotyping/$region/CoordenadasRelativas/C
              > $wd/Genotyping/$region/Sondas/C
          elif [ ! -f $wd/Genotyping/$region/Sondas/D ];then
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
    fi
done

#rm -r $wd/Infor/probecreation
#rm -r $wd/Infor/windows
#rm -r $wd/Infor/blast_results
