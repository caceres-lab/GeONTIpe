#!/bin/bash

wd=$(pwd)

Ind="$(cut -f1 $wd/tracking_pipeline/completed_samples.txt)"
Invs="$(cut -f1 $wd/Infor/ListaRef.txt)"
Lista="$wd/Infor/ListaRef.txt"

for i in $Ind
do

  echo -e "Inv\tGenotype" > $wd/${i}_Genotype.txt
  echo -e "Inv\tSV" > $wd/${i}_SV.txt
  echo -e "Inv\tReads" > $wd/${i}_reads.txt

 for inv in $Invs
 do
    if [ -f $wd/${i}/Genotyping/${inv}/${inv}_Genotype.txt ];then
      countStd="$(cat $wd/${i}/Genotyping/${inv}/${inv}_Genotype.txt | tail -n +2 | cut -f3 | grep -c "Std")"
      countInv="$(cat $wd/${i}/Genotyping/${inv}/${inv}_Genotype.txt | tail -n +2 | cut -f3 | grep -c "Inv")"
      count="$((countStd+countInv))"
        if [ -d $wd/${i}/Genotyping/${inv}/SNP ];then
            if [ -f $wd/${i}/Genotyping/${inv}/SNP/${inv}_resolvedA.txt ] || [ -f $wd/${i}/Genotyping/${inv}/SNP/${inv}_resolvedD.txt ];then
                if [ "$(cat $wd/${i}/Genotyping/${inv}/${inv}_Genotype.txt | tail -n +2 | cut -f3 | grep -c "Std")" -ge 1 ];then
                  Geno="Std"
                elif [ "$(cat $wd/${i}/Genotyping/${inv}/${inv}_Genotype.txt | tail -n +2 | cut -f3 | grep -c "Inv")" -ge 1 ];then
                  Geno="Inv"
                fi
            else
              Geno="NER"
            fi                
        else
          Geno="$(cat $wd/${i}/Genotyping/${inv}/${inv}_Genotype.txt | grep FinalGenotype | cut -f2)"
        fi
    else
      Geno="ND"
      count=0
    fi

    if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_AllBPs.txt ];then
          SV1=""
          SV2=""
          SV3=""
          SV4=""
          SV5=""
          SV6=""
          SV1sx1=""
          SV1sx2=""
          SV2sx1=""
          SV2sx2=""
          SV3sx1=""
          SV3sx2=""
          SV4sx1=""
          SV4sx2=""
          SV5sx1=""
          SV5sx2=""
          SV6sx1=""
          SV6sx2=""
  
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP1_Std.txt ];then
            SVs1="$(cat $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP1_Std.txt | tail -n +2 | cut -f2,3 | sort | uniq | sed 's/\t/:/g')"
            bpsv1="$(ls $wd/${i}/Genotyping/${inv} | grep SVBP1_Std | sed 's/_/\t/'g | cut -f3 | sed 's/SV//g' | sed 's/BP1/BP1;/g')"		
                    
              if [ "$(echo $SVs1 | sed 's/ /\n/g' | wc -l)" == 2 ];then
                      SV1sx1="$(echo $SVs1 | sed 's/ /\n/g' | head -1)"
                SV1sx2="$(echo $SVs1 | sed 's/ /\n/g' | head -2 | tail -1)"
                SV1="${SV1sx1}"_"${bpsv1}""${SV1sx2}"_"${bpsv1}"
                    elif [ "$(echo $SVs1 | sed 's/ /\n/g' | wc -l)" -ge 3 ];then
                SV1="error"
              else
                SV1="${SVs1}"_"${bpsv1}"
              fi
          fi
          
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP1_Inv.txt ];then
            SVs2="$(cat $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP1_Inv.txt | tail -n +2 | cut -f2,3 | sort | uniq | sed 's/\t/:/g')"
            bpsv2="$(ls $wd/${i}/Genotyping/${inv} | grep SVBP1_Inv | sed 's/_/\t/'g | cut -f3 | sed 's/SV//g' | sed 's/BP1/BP1;/g')"
            
            if [ "$(echo $SVs2 | sed 's/ /\n/g' | wc -l)" == 2 ];then
                      SV2sx1="$(echo $SVs2 | sed 's/ /\n/g' | head -1 )"
                SV2sx2="$(echo $SVs2 | sed 's/ /\n/g' | head -2 | tail -1)"
                SV2="${SV2sx1}"_"${bpsv2}""${SV2sx2}"_"${bpsv2}"
                  elif [ "$(echo $SVs2 | sed 's/ /\n/g' | wc -l)" -ge 3 ];then
                SV2="error"
            else
                SV2="${SVs2}"_"${bpsv2}"
              fi
          fi
          
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP2_Std.txt ];then
            SVs3="$(cat $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP2_Std.txt | tail -n +2 | cut -f2,3 | sort | uniq | sed 's/\t/:/g')"
            bpsv3="$(ls $wd/${i}/Genotyping/${inv} | grep SVBP2_Std | sed 's/_/\t/'g | cut -f3 | sed 's/SV//g' | sed 's/BP2/BP2;/g')"
            
            if [ "$(echo $SVs3 | sed 's/ /\n/g' | wc -l)" == 2 ];then
                      SV3sx1="$(echo $SVs3 | sed 's/ /\n/g' | head -1)"
                SV3sx2="$(echo $SVs3 | sed 's/ /\n/g' | head -2 | tail -1)"
                SV3="${SV3sx1}"_"${bpsv3}""${SV3sx2}"_"${bpsv3}"
                  elif [ "$(echo $SVs3 | sed 's/ /\n/g' | wc -l)" -ge 3 ];then
                SV3="error"
            else
                SV3="${SVs3}"_"${bpsv3}"
            fi
          fi
          
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP2_Inv.txt ];then
            SVs4="$(cat $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBP2_Inv.txt | tail -n +2 | cut -f2,3 | sort | uniq | sed 's/\t/:/g')"
            bpsv4="$(ls $wd/${i}/Genotyping/${inv} | grep SVBP2_Inv | sed 's/_/\t/'g | cut -f3 | sed 's/SV//g' | sed 's/BP2/BP2;/g')"
            
            if [ "$(echo $SVs4 | sed 's/ /\n/g' | wc -l)" == 2 ];then
                      SV4sx1="$(echo $SVs4 | sed 's/ /\n/g' | head -1 )"
                SV4sx2="$(echo $SVs4 | sed 's/ /\n/g' | head -2 | tail -1)"
                SV4="${SV4sx1}"_"${bpsv4}""${SV4sx2}"_"${bpsv4}"
                  elif [ "$(echo $SVs4 | sed 's/ /\n/g' | wc -l)" -ge 3 ];then
                SV4="error"
            else
                SV4="${SVs4}"_"${bpsv4}"
            fi
          fi
          
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBC.txt ];then
            SVs5="$(cat $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVBC.txt | tail -n +2 | cut -f2,3 | sort | uniq | sed 's/\t/:/g')"
            bpsv5="$(ls $wd/${i}/Genotyping/${inv} | grep SVBC | sed 's/_/\t/'g | cut -f3 | sed 's/SV//g' | sed 's/BC.txt/BC;/g')"
            
            if [ "$(echo $SVs5 | sed 's/ /\n/g' | wc -l)" == 2 ];then
                      SV5sx1="$(echo $SVs5 | sed 's/ /\n/g' | head -1)"
                SV5sx2="$(echo $SVs5 | sed 's/ /\n/g' | head -2 | tail -1)"
                SV5="${SV5sx1}"_"${bpsv5}""${SV5sx2}"_"${bpsv5}"
                  elif [ "$(echo $SVs5 | sed 's/ /\n/g' | wc -l)" -ge 3 ];then
                SV5="error"
            else
                SV5="${SVs5}"_"${bpsv5}"
              fi
          fi
          
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVAD.txt ];then
            SVs6="$(cat $wd/${i}/Genotyping/${inv}/DetSV_${inv}_SVAD.txt | tail -n +2 | cut -f2,3 | sort | uniq | sed 's/\t/:/g')"
            bpsv6="$(ls $wd/${i}/Genotyping/${inv} | grep SVAD | sed 's/_/\t/'g | cut -f3 | sed 's/SV//g' | sed 's/AD.txt/AD;/g')"
            
            if [ "$(echo $SVs6 | sed 's/ /\n/g' | wc -l)" == 2 ];then
                      SV6sx1="$(echo $SVs6 | sed 's/ /\n/g' | head -1 )"
                SV6sx2="$(echo $SVs6 | sed 's/ /\n/g' | head -2 | tail -1)"
                SV6="${SV6sx1}"_"${bpsv6}""${SV6sx2}"_"${bpsv6}"
                  elif [ "$(echo $SVs6 | sed 's/ /\n/g' | wc -l)" -ge 3 ];then
                SV6="error"
            else
                SV6="${SVs6}"_"${bpsv6}"
            fi
          fi
          
          if [ -f $wd/${i}/Genotyping/${inv}/DetSV_${inv}_AllBPs.txt ];then
            SVs="${SV1}${SV2}${SV3}${SV4}${SV5}${SV6}"
          else
            SVs="-"
          fi
    else
        SVs="-"
    fi
  
    echo -e "$inv\t$Geno" >> $wd/${i}_Genotype.txt
    echo -e "$inv\t$SVs" >> $wd/${i}_SV.txt
    echo -e "$inv\t$count" >> $wd/${i}_reads.txt
 done
done

mkdir -p $wd/Genos
mkdir -p $wd/SVs
mkdir -p $wd/counts

mv $wd/*_Genotype.txt $wd/Genos
mv $wd/*_SV.txt $wd/SVs
mv $wd/*_reads.txt $wd/counts

echo "invs" > "$wd/Genotypes.tmp"
cut -f1 "$Lista" >> "$wd/Genotypes.tmp"

echo "invs" > "$wd/SV.tmp"
cut -f1 "$Lista" >> "$wd/SV.tmp"

echo "invs" > "$wd/count.tmp"
cut -f1 ""$Lista"" >> "$wd/count.tmp"

for a in $Ind 
do
  individual="$a"
  genotypes="$(cut -f2 "$wd/Genos/${a}_Genotype.txt" | tail -n +2)"
  cat <(echo -e "$individual") <(echo -e "$genotypes") > "$wd/allgenos_$a.tmp"

  paste -d '\t' "$wd/Genotypes.tmp" "$wd/allgenos_$a.tmp" > "$wd/Genotypes_tmp.tmp"
  mv "$wd/Genotypes_tmp.tmp" "$wd/Genotypes.tmp"
  
  svs="$(cut -f2 "$wd/SVs/${a}_SV.txt" | tail -n +2)"
  cat <(echo -e "$individual") <(echo -e "$svs") > "$wd/allsvs_$a.tmp"  
  
  paste -d '\t' "$wd/SV.tmp" "$wd/allsvs_$a.tmp" > "$wd/svs_tmp.tmp"
  mv "$wd/svs_tmp.tmp" "$wd/SV.tmp"
  
  reads="$(cut -f2 "$wd/counts/${a}_reads.txt" | tail -n +2)"
  cat <(echo -e "$individual") <(echo -e "$reads") > "$wd/allreads_$a.tmp"
  
  paste -d '\t' "$wd/count.tmp" "$wd/allreads_$a.tmp" > "$wd/reads_tmp.tmp"
  mv "$wd/reads_tmp.tmp" "$wd/count.tmp"
done

mv $wd/Genotypes.tmp $wd/Genotypes.txt
mv $wd/SV.tmp $wd/SV.txt
mv $wd/count.tmp $wd/Reads.txt

rm "$wd/allgenos"*.tmp
rm "$wd/allsvs"*.tmp
rm "$wd/allreads"*.tmp



