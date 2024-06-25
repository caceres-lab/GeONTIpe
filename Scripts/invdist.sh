#!/bin/bash

wd="$(pwd)"

echo -e "Inv\tDistAC\tDistBD\tDistAB\tDistCD" > $wd/Infor/ExpectedDistProvesInverted.txt

for i in $(cut -f1 $wd/Infor/ListaRef.txt)
do

#### Calculate the mean position of the probe

echo $i

	A1="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/A | cut -f1)"
	A2="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/A | cut -f2 | tr -d '\r')"
	B1="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/B | cut -f1)"
	B2="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/B | cut -f2 | tr -d '\r')"
	C1="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/C | cut -f1)"
	C2="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/C | cut -f2 | tr -d '\r')"
	D1="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/D | cut -f1)"
	D2="$(cat $wd/Genotyping/${i}/CoordenadasRelativas/D | cut -f2 | tr -d '\r')"

	sumA=$((A1 + A2))
	sumB=$((B1 + B2))
	sumC=$((C1 + C2))
	sumD=$((D1 + D2))

	A=$(($sumA/2))
	B=$(($sumB/2))
	C=$(($sumC/2))
	D=$(($sumD/2))

#### Add the absolute position

	Inicio="$(grep ${i}[[:space:]] $wd/Infor/ListaRef.txt | sed 's/:/\t/g' | sed 's/-/\t/g' |  cut -f3 )"
 
	Ai=$((A + Inicio))
	Bi=$((B + Inicio))
	Ci=$((C + Inicio))
	Di=$((D + Inicio))
 
	BP1A="$(grep ${i}[[:space:]] $wd/Infor/allcoords.txt | cut -f3)"
	BP1B="$(grep ${i}[[:space:]] $wd/Infor/allcoords.txt | cut -f4)"
	BP2A="$(grep ${i}[[:space:]] $wd/Infor/allcoords.txt | cut -f5)"
	BP2B="$(grep ${i}[[:space:]] $wd/Infor/allcoords.txt | cut -f6 | tr -d '\r')"

#### Calculate the distances between probes

	DistA=$(echo "$BP1B-$Ai" | bc) 
	DistA_abs=${DistA#-}
	DistC1=$(echo "$BP2A-$Ci" | bc) 
	DistC1_abs=${DistC1#-}
	DistC2=$(echo "$Ci-$BP1B" | bc)
	DistC2_abs=${DistC2#-}
	DistD=$(echo "$Di-$BP2A" | bc)
	DistD_abs=${DistD#-}
	DistB1=$(echo "$Bi-$BP1B" | bc)
	DistB1_abs=${DistB1#-}
	DistB2=$(echo "$BP2A-$Bi" | bc) 
	DistB2_abs=${DistB2#-}

	DistAC=$((DistA_abs + DistC1_abs))
	DistBD=$((DistD_abs + DistB1_abs))
	DistAB=$((DistA_abs + DistB2_abs))
	DistCD=$((DistD_abs + DistC2_abs))

	echo -e "${i}\t$DistAC\t$DistBD\t$DistAB\t$DistCD" >> $wd/Infor/ExpectedDistProvesInverted.txt
done
