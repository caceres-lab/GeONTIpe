#!/bin/bash

## Parameters
part="$1"
quality="$2"
length="$3"
reads="$4"
inv_det="$5"
sv_det="$6"
reference="$7"
mapq="$8"
wd="$9"

sp="$(echo $part | sed 's/.part/\t/g' | cut -f1)"
split="$(echo $sp.fastq.gz)"
individuo="$(echo $part | sed 's/.part/\t/g' | sed 's/_/\t/g' | cut -f1)"

## Filtering

mkdir -p $wd/$individuo/Resultados/mappeo

if [ -f $wd/$individuo/Resultados/mappeo/Sorted${part}.bam ] && [ -s $wd/$individuo/Resultados/mappeo/Sorted${part}.bam ];then
  continue
elif [ ! -f $wd/$individuo/Resultados/Filt_fq/Filt_${part}.fastq.gz ] || [ ! -s $wd/$individuo/Resultados/Filt_fq/Filt_${part}.fastq.gz ];then
  gunzip -c $wd/$individuo/Descargas/$individuo/${split}.split/${part}.fastq.gz | NanoFilt  -l $length -q $quality | gzip > $wd/$individuo/Resultados/Filt_fq/Filt_${part}.fastq.gz 
fi

## Mapping

if [ -f $wd/$individuo/Resultados/mappeo/Sorted${part}.bam ] && [ -s $wd/$individuo/Resultados/mappeo/Sorted${part}.bam ];then
    continue
    echo -e "$individuo\t$part\tcompleted" >> $wd/tracking_pipeline/progress.txt
else
  minimap2 -ax $reads -z $inv_det -r $sv_det --secondary=no $reference $wd/$individuo/Resultados/Filt_fq/Filt_${part}.fastq.gz > $wd/$individuo/Resultados/mappeo/Map_${part}.bam 
  samtools view -bS $wd/$individuo/Resultados/mappeo/Map_${part}.bam | samtools view -b -F 0x4 - --input-fmt-option "filter= mapq>=${mapq}" > $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam
  samtools view -bS $wd/$individuo/Resultados/mappeo/Map_${part}.bam | samtools view -b -f 0x4 - > $wd/$individuo/Resultados/mappeo/Unmapped_${part}.bam
  rm $wd/$individuo/Resultados/mappeo/Map_${part}.bam
  if [ -f $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam ] && [ -s $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam ];then
    samtools sort -o $wd/$individuo/Resultados/mappeo/Sorted_${part}.bam $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam
    samtools sort -o $wd/$individuo/Resultados/mappeo/Unmappedsorted_${part}.bam $wd/$individuo/Resultados/mappeo/Unmapped_${part}.bam
    rm $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam
    rm $wd/$individuo/Resultados/mappeo/Unmapped_${part}.bam
    rm $wd/$individuo/Resultados/Filt_fq/Filt_${part}.fastq.gz
    rm $wd/$individuo/Descargas/$individuo/${split}.split/${part}.fastq.gz
    echo -e "$individuo\t$part\tcompleted" >> $wd/tracking_pipeline/progress.txt
  elif [ ! -f $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam ] || [ ! -s $wd/$individuo/Resultados/mappeo/Mapped_${part}.bam ];then
    echo -e "$individuo\t$part\tfailed" >> $wd/tracking_pipeline/progress.txt
  fi
fi
