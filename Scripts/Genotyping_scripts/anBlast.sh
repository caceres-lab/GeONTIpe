#!/bin/bash

wd="$1"

inv="$2"
task="$3"
cov_probe="$4"
id_probe="$5"

if [ $# -eq 0 ]
then
	echo "No arguments"
	exit
fi

mkdir -p $wd/$inv/Resultados/Distancias

samtools fastq $wd/$inv/${inv}DupReads.bam > $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq

> $wd/$inv/Resultados/sondaErronea

for dread in $(cat $wd/$inv/Resultados/readsDup)
do
	mkdir -p $wd/$inv/Resultados/Distancias/Reads/$dread

	numfila=$(cat -n $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq | grep $dread | cut -f1 | tr -d " ")
	tail -n+$numfila $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq | head -4 | head -1 | sed 's/@/#/g' > $wd/$inv/Resultados/Distancias/Reads/$dread/nameread.fastq
  tail -n+$numfila $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq | head -4 | tail -3 > $wd/$inv/Resultados/Distancias/Reads/$dread/seqread.fastq
  cat $wd/$inv/Resultados/Distancias/Reads/$dread/nameread.fastq $wd/$inv/Resultados/Distancias/Reads/$dread/seqread.fastq > $wd/$inv/Resultados/Distancias/Reads/$dread/read.fastq
  rm $wd/$inv/Resultados/Distancias/Reads/$dread/nameread.fastq
  rm $wd/$inv/Resultados/Distancias/Reads/$dread/seqread.fastq

	## Loop for every probe

	for sonda in A B C D
	do
		mkdir -p $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord
		blastn -query $wd/$inv/Sondas/${sonda} -subject $wd/$inv/Resultados/Distancias/Reads/${dread}/read.fastq -task $task -outfmt "6 sstart send slen sstrand" -qcov_hsp_perc $cov_probe -perc_identity $id_probe -sorthits 4 > $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda}

		if [ $(cat $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda} | wc -l) -eq 0 ]
		then
			rm $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda}
		fi

		if [ $(cat $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda} | wc -l) -gt 1 ]
		then
			echo -e "${dread}\t${sonda}" >> $wd/$inv/Resultados/sondaErronea
		fi

	done

done
