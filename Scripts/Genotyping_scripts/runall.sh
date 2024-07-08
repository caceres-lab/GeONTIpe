#!/bin/bash

wd="$1"

invs="$(cut -f1 $wd/../../Infor/ListaRef.txt)"

ind=$2
local=$3
route=$4

samplgender="$(cat $wd/../../Infor/gender | grep $ind[[:space:]] | cut -f2 | tr -d '\r')"

if [ $samplgender == "Women" ];then
  gender="W"
elif [ $samplgender == "Men" ];then
  gender="M"
fi

echo ${gender}

if [ $local == Y ];then
  if [ -f "$wd/../../../cosas/"$ind"/hg38/"$ind".bam" ] && [ -f "$wd/../Descargas/Resultados/mappeo/"$ind".bam" ]; then
    initBam="$wd/../../../cosas/"$ind"/hg38/"$ind".bam"
    initBam2="$wd/../Descargas/Resultados/mappeo/"$ind".bam"
  elif [ -f "$wd/../../../cosas/$ind/hg38/"$ind".bam" ]; then
    initBam="$wd/../../../cosas/"$ind"/hg38/"$ind".bam"
  else
    initBam="$wd/../Resultados/mappeo/"$ind".bam"
  fi
else
  initBam=$route/$ind".hg38.cram"
fi

for inv in $invs
do

  if [ -f $wd/$inv/${inv}.png ] && [ -d $wd/$inv/Resultados ];then
		continue
  fi
  
  if [ ! -f $wd/$inv/${inv}.png ] && [ -d $wd/$inv/Resultados ] && [ -s $wd/$inv/Resultados/MapInfo ]; then
    Rscript $wd/signos.R ${wd} ${inv} ${gender} ${8} ${9} ${10} ${11} ${12}
    continue
	fi

	echo $inv
 
 ## Structure of the inversion should have the following chr:pos1-pos2, if is not encoded in that way you should add extra steps to you input be coherent with your mapping files

 chr=$(grep "$inv[[:space:]]" $wd/../../Infor/ListaRef.txt | cut -f2 | sed 's/:/\t/g' | cut -f1)
 #chr_name="$(cat $wd/../../Infor/conversion.txt | grep $chr[[:space:]] | cut -f2)"
 #chr_name=chr$chr
 pos=$(grep "$inv[[:space:]]" $wd/../../Infor/ListaRef.txt | cut -f2 | tr -d "[[:space:]]" | sed 's/:/\t/g' | cut -f2 | tr -d "[[:space:]]")
 #selInv=$chr_name:$pos
 selInv=$chr:$pos
 echo $selInv | sed 's/chr//; s/:/\t/g; s/-/\t/g' | sed 's/Y/24/g' | sed 's/X/23/g' > $wd/$inv/seqInfo
  
  ## Extraction of the region from the main file

	if [ -f "$wd/../../../cosas/"$ind"/hg38/FinalMap.bam" ] && [ -f "$wd/../Descargas/Resultados/mappeo/FinalMap.bam" ];then
		samtools view -b $initBam "$selInv" > $wd/$inv/${inv}_1.bam
		samtools view -b $initBam2 "$selInv" > $wd/$inv/${inv}_2.bam
		samtools merge $wd/$inv/${inv}.bam $wd/$inv/${inv}_1.bam $wd/$inv/${inv}_2.bam
		samtools index $wd/$inv/$inv.bam > $wd/$inv/$inv.bam.bai
	else
		samtools view -b $initBam "$selInv" > $wd/$inv/${inv}.bam
		samtools index $wd/$inv/$inv.bam > $wd/$inv/$inv.bam.bai
	fi
  
  if [ $local == N ];then
    rm $wd/../../*.cram.crai
	fi
  
  ## Mapping the reads agains the probes

	mkdir -p $wd/$inv/Resultados/SondasMapeadas
	letras="A B C D"

	> $wd/$inv/Resultados/MapInfo
 
   laInit="$(($(cut -f2 $wd/$inv/seqInfo)-1))"

	for sonda in $letras
	do
		initSonda="$(($(cut -f1 $wd/$inv/CoordenadasRelativas/$sonda | tr -d "[[:space:]]")+$laInit))"
		finSonda="$(($(cut -f2 $wd/$inv/CoordenadasRelativas/$sonda | tr -d "[[:space:]]")+$laInit))"

		samtools view $wd/$inv/$inv.bam -F 0x0010 --input-fmt-option "filter= pos <= $initSonda && endpos >= $finSonda" | cut -f1 > $wd/$inv/Resultados/SondasMapeadas/${sonda}_f
		samtools view $wd/$inv/$inv.bam -f 0x0010 --input-fmt-option "filter= pos <= $initSonda && endpos >= $finSonda" | cut -f1 > $wd/$inv/Resultados/SondasMapeadas/${sonda}_r

		> $wd/$inv/Resultados/SondasMapeadas/${sonda}

		for line in $(cat $wd/$inv/Resultados/SondasMapeadas/${sonda}_f)
		do
			echo -e "$line\t+" >> $wd/$inv/Resultados/SondasMapeadas/${sonda}
			echo -e "$line\t$sonda\t+" >> $wd/$inv/Resultados/MapInfo
		done

		for line in $(cat $wd/$inv/Resultados/SondasMapeadas/${sonda}_r)
		do
			echo -e "$line\t-" >> $wd/$inv/Resultados/SondasMapeadas/${sonda}
			echo -e "$line\t$sonda\t-" >> $wd/$inv/Resultados/MapInfo
		done

		rm $wd/$inv/Resultados/SondasMapeadas/${sonda}_*

	done

  ## Generating a table with the hits of the relative coordinates of the probes
	todosReads=$(cat $wd/$inv/Resultados/SondasMapeadas/* | cut -f1 | sort | uniq)
	> $wd/$inv/Resultados/TablaReadTF

	for mapRead in $todosReads
	do
		siA="F"
		siB="F"
		siC="F"
		siD="F"

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/A ;
		then
			siA="T"
		fi

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/B ;
		then
			siB="T"
		fi

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/C ;
		then
			siC="T"
		fi

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/D ;
		then
			siD="T"
		fi

		echo -e "$mapRead\t$siA\t$siB\t$siC\t$siD" >> $wd/$inv/Resultados/TablaReadTF

	done

	for sonda in $letras; do cut -f1 $wd/$inv/Resultados/SondasMapeadas/$sonda | sort | uniq; done | sort | uniq > $wd/$inv/Resultados/readsDup

	cat $wd/$inv/Resultados/TablaReadTF | grep -f $wd/$inv/Resultados/readsDup > $wd/$inv/Resultados/TablaCandidatos

	cut -f1 $wd/$inv/Resultados/TablaCandidatos > $wd/$inv/Resultados/readsDup

	nUnaSonda="$(cut -f1 $wd/$inv/Resultados/MapInfo | sort | uniq | wc -l)"
	nDupSonda="$(cut -f1 $wd/$inv/Resultados/MapInfo | sort | uniq -d | wc -l)"

	if [ $nUnaSonda -gt 300 ] && [ $nDupSonda -gt 10 ]
	then
		cut -f1 $wd/$inv/Resultados/MapInfo | sort | uniq -d > $wd/$inv/Resultados/readsDup
		cat $wd/$inv/Resultados/TablaReadTF | grep -f $wd/$inv/Resultados/readsDup | awk '{if ($2 == "F" && $3 == "T" && $4 == "T" && $5 == "F"); else print}' | cut -f1 > $wd/$inv/Resultados/tmp
		mv $wd/$inv/Resultados/tmp $wd/$inv/Resultados/readsDup
	fi

	samtools view -H $wd/$inv/$inv.bam > $wd/$inv/header
	samtools view $wd/$inv/$inv.bam | grep -f $wd/$inv/Resultados/readsDup > $wd/$inv/reads

	cat $wd/$inv/header $wd/$inv/reads | samtools view -bS -o $wd/$inv/${inv}DupReads.bam
	rm $wd/$inv/header $wd/$inv/reads

	echo "Analyzing BLAST of $inv"
	$wd/anBlast.sh ${wd} ${inv} ${5} ${6} ${7}

  Rscript $wd/signos.R $wd ${inv} ${gender} ${8} ${9} ${10} ${11} ${12}
  
  if [ -f $wd/$inv/temp_SNP.txt ]; then
    echo "Analyzing SNPs of $inv"
    #rm "$wd/$inv/temp_SNP.txt"
    $wd/snps.sh ${inv} ${2} ${13} ${14} ${15} ${16} ${17} ${18}
    rm $wd/$inv/Regions
  fi

	rm $wd/$inv/${inv}.bam
	rm $wd/$inv/${inv}.bam.bai
	samtools index $wd/$inv/${inv}DupReads.bam > $wd/$inv/${inv}DupReads.bam.bai
  rm $wd/$inv/${inv}DupReads.bam
  rm $wd/$inv/${inv}DupReads.bam.bai
  rm $wd/$inv/Regions
done
