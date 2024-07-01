#!/bin/bash

wd="$(pwd)"

> $wd/Infor/probecreation/allregions.txt

for region in $(cut -f1 Infor/allcoords.txt)
do
    if [ -d $wd/Genotyping/$region ] && [ $(ls $wd/Genotyping/$region/Sondas | wc -l) == 4 ]; then
        continue
    else
        if [ -d $wd/Infor/windows/$region/A ]; then
            ls $wd/Infor/windows/$region/A | grep -v "allwindows" | sort -rn > $wd/Infor/windows/$region/A/allwindows.txt
            echo $region/A >> $wd/Infor/probecreation/allregions.txt
        fi

        if [ -d $wd/Infor/windows/$region/B ]; then
            ls $wd/Infor/windows/$region/B | grep -v "allwindows" | sort -n > $wd/Infor/windows/$region/B/allwindows.txt
            echo $region/B >> $wd/Infor/probecreation/allregions.txt
        fi

        if [ -d $wd/Infor/windows/$region/C ]; then
            ls $wd/Infor/windows/$region/C | grep -v "allwindows" | sort -rn > $wd/Infor/windows/$region/C/allwindows.txt
            echo $region/C >> $wd/Infor/probecreation/allregions.txt
        fi

        if [ -d $wd/Infor/windows/$region/D ]; then
            ls $wd/Infor/windows/$region/D | grep -v "allwindows" | sort -n > $wd/Infor/windows/$region/D/allwindows.txt
            echo $region/D >> $wd/Infor/probecreation/allregions.txt
        fi
    fi
done
