#!/bin/bash

inv=$1
ind=$2
elsnp=$3
ref=$4
quality=$5

wd="$(pwd)/$ind/Genotyping"

initBam="$wd/$inv/${inv}.bam"

> $wd/Snps_review.txt
> $wd/snp_resolution.txt

if [ -f $wd/$inv/SNP/SNPs_${inv}.png ]; then 
  continue
elif [ -f $wd/$inv/SNP/variante${inv} ]; then
  if [ "$(cat $wd/$inv/SNP/variante${inv} | wc -l)" != "0" ]; then
    Rscript $wd/parsempu.R $wd/$inv/SNP/variante${inv} $wd/$inv/${inv}_Genotype.txt $wd/$inv/SNP/snps${inv}.txt  $wd/$inv/Regions ${inv} ${6} ${7} ${8} $wd
  fi
else
  
  mkdir -p $wd/$inv/SNP

# Selection of SNPs in the positions

  samtools faidx $ref
  
  chr=$(cut -f1 $wd/$inv/seqInfo)
  if [ $chr == "23" ]
  then
    chr="X"
  elif [ $chr == "24" ]
  then
    chr="Y"
  fi
  
  > $wd/$inv/Regions
  
  A="$(grep "${inv}[[:space:]]" $wd/../../Infor/allcoords.txt | cut -f3)"
  B="$(grep "${inv}[[:space:]]" $wd/../../Infor/allcoords.txt | cut -f4)"
  C="$(grep "${inv}[[:space:]]" $wd/../../Infor/allcoords.txt | cut -f5)"
  D="$(grep "${inv}[[:space:]]" $wd/../../Infor/allcoords.txt | cut -f6)"
  pos1="$(grep "${inv}[[:space:]]" $wd/../../Infor/ListaRef.txt | sed 's/:/\t/g' | cut -f3 | sed 's/-/\t/g' | cut -f1)"
  pos2="$(grep "${inv}[[:space:]]" $wd/../../Infor/ListaRef.txt | sed 's/:/\t/g' | cut -f3 | sed 's/-/\t/g' | cut -f2 | tr -d '\r' )"
  
  echo -e "${inv}_A\t$chr\t$pos1\t$A" >> $wd/$inv/Regions
  echo -e "${inv}_BC\t$chr\t$B\t$C" >> $wd/$inv/Regions
  echo -e "${inv}_D\t$chr\t$D\t$pos2" >> $wd/$inv/Regions
  
  bcftools view -r ${chr}:$pos1-$A -m2 -M2 -v snps ${elsnp} | grep -v "#" > $wd/$inv/SNP/snps${inv}.txt
  bcftools view -r ${chr}:$B-$C -m2 -M2 -v snps ${elsnp} | grep -v "#" >> $wd/$inv/SNP/snps${inv}.txt
  bcftools view -r ${chr}:$D-$pos2 -m2 -M2 -v snps ${elsnp} | grep -v "#" >> $wd/$inv/SNP/snps${inv}.txt

### See the bases of X quality on the regions selected in map file
   
  samtools mpileup --no-BAQ -f $ref -Q $quality --output-QNAME $initBam > $wd/$inv/SNP/mpileup${inv}.txt
  
#### Check if sample exists in 1000 Genome project

  chrs="$(echo chr${chr})"

if [ $ind == "HG001" ];then
  ind="NA12878"
elif [ $ind == "HG002" ];then
  ind="NA24385"
elif [ $ind == "HG003" ];then
  ind="NA24149"
elif [ $ind == "HG004" ];then
  ind="NA24143"
elif [ $ind == "HG005" ];then
  ind="NA24631"
elif [ $ind == "HG006" ];then
  ind="NA24694"
elif [ $ind == "HG007" ];then
  ind="NA24695"  
fi

if [[ $(bcftools view -h "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${chrs}.recalibrated_variants.vcf.gz" | grep -c "${ind}") == 1 ]]; then
  bcftools view -H -r ${chrs}:$pos1-$pos2 -s "${ind}" http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${chrs}.recalibrated_variants.vcf.gz | cut -f2,10 | sed 's/:/\t/g' | cut -f1-2 | grep "0/1" | cut -f1 > $wd/$inv/SNP/snps1000GP${inv}.txt
fi

### Search the positions in heterozygosis in 1000 Genome Project and filter out those with no heterozygosis   

if [ -f $wd/$inv/SNP/snps1000GP${inv}.txt ];then
  mv $wd/$inv/SNP/snps${inv}.txt $wd/$inv/SNP/snpslack1000GP${inv}.txt
  
  grep -w -f $wd/$inv/SNP/snps1000GP${inv}.txt $wd/$inv/SNP/snpslack1000GP${inv}.txt > $wd/$inv/SNP/snps${inv}.txt
fi

### Select positions of mpileup with selected variants

>$wd/$inv/SNP/variante${inv}
  
  for snp in $(cat $wd/$inv/SNP/snps${inv}.txt | cut -f2)
  do
    cat $wd/$inv/SNP/mpileup${inv}.txt | grep "[[:space:]]$snp[[:space:]]" | cut -f2,3,5,7 >> $wd/$inv/SNP/variante${inv}
  done

#rm $wd/$inv/SNP/mpileup${inv}.txt

if [ ! -f $wd/$inv/SNP/variante${inv} ] && [ ! -f $wd/$inv/SNP/snps${inv}.txt ]; then
  echo "Not variants found on the postions" >> $wd/$inv/errors.txt
elif [ ! -s $wd/$inv/SNP/variante${inv} ] && [ ! -s $wd/$inv/SNP/snps${inv}.txt ]; then
  echo "Not variants found on the postions" >> $wd/$inv/errors.txt
else
  Rscript $wd/parsempu.R $wd/$inv/SNP/variante${inv} $wd/$inv/${inv}_Genotype.txt $wd/$inv/SNP/snps${inv}.txt $wd/$inv/Regions ${inv} ${6} ${7} ${8} $wd
fi

fi

if [ "$(ls $wd/${inv}/SNP | grep -c "error")" != 0 ];then
  cat $wd/${inv}/SNP/${inv}_error*.txt > $wd/${inv}/SNP/${inv}_fails.txt
  rm $wd/${inv}/SNP/${inv}_error*.txt
fi
  
if [ -f $wd/$inv/SNP/${inv}_resolvedA.txt ] && [ -f $wd/$inv/SNP/${inv}_resolvedD.txt ];then
  echo "${inv} resolved both" >> $wd/snp_resolution.txt
elif [ -f $wd/$inv/SNP/${inv}_resolvedA.txt ]; then
  echo "${inv} resolved A" >> $wd/snp_resolution.txt
elif [ -f $wd/$inv/SNP/${inv}_resolvedD.txt ]; then
  echo "${inv} resolved D" >> $wd/snp_resolution.txt
elif [ -f $wd/$inv/SNP/${inv}_not_resolved_A.txt ] && [ -f $wd/$inv/SNP/${inv}_not_resolvedD.txt ];then
  echo "${inv} not resolt" >> $wd/snp_resolution.txt
else
  echo "${inv} not resolt" >> $wd/snp_resolution.txt
fi

rm "$(pwd)"/*.vcf.gz.tbi

