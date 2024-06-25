#!/bin/bash

wd="$(pwd)"

inv="$1"

if [ -f $wd/$inv/Sondas/A.txt ]
then
	for sonda in A B C D
	do
		mv $wd/$inv/Sondas/${sonda}.txt $wd/$inv/Sondas/${sonda}
		mv $wd/$inv/CoordenadasRelativas/${sonda}.txt $wd/$inv/CoordenadasRelativas/${sonda}
	done

fi
