#!/bin/bash

wd="$(pwd)"

## Parametros configurables
ind="$1"
separate="$2"

files="$(ls $wd/${ind}/Descargas/${ind} | sed 's/.fastq.gz//g')"

for split in $files
do
  seqkit split2 -p $separate $wd/${ind}/Descargas/${ind}/$split.fastq.gz 
  ls $wd/${ind}/Descargas/${ind}/${split}.fastq.gz.split >> $wd/${ind}/Descargas/${ind}/List_${individuo}
  cat $wd/${ind}/Descargas/${ind}/List_${individuo} | sed 's/.fastq.gz//g' >> $wd/FileList.txt
  rm $wd/${ind}/Descargas/${ind}/List_${individuo}
done
