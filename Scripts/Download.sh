#!/bin/bash

wd="$(pwd)"

## Parametros configurables
ind="$1"
fastq="$(grep ${ind} $wd/Infor/ListaTodos.txt | cut -f2 | tr -d '\r')"

contador=1
for secuencia in $fastq
do
   name_file="${ind}_${contador}.fastq.gz" 
   wget -P $wd/${ind}/Descargas/${ind} -O $wd/${ind}/Descargas/${ind}/$name_file $secuencia
   contador=$((contador + 1))
done

ls $wd/${ind}/Descargas/${ind} | sed 's/.fastq.gz//g' >> $wd/tracking_pipeline/Download_fastq.txt