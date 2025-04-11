#!/bin/bash

problem=$1
InitialGuesses=$2
nvert=$3

if [ -z "$problem" ] || [ -z "$InitialGuesses" ] || \
   { [[ "$problem" == "regpol" ]] && [[ -z "$nvert" ]];}; then
    echo "Usage: $0 Agen|regpol random|lattice [nvert]"
    exit 1
fi

make

backupdir=~/documents/covering/tests/MIF${InitialGuesses}_$(date +%Y%m%d)
balls=$(seq 1 1 10)
#balls=8
budget=(0 100 200 500 1000)

if [ ! -d "$backupdir" ]; then #Check whether the variable is an empty string
    mkdir -p $backupdir
fi

for ((i=1; i<${#budget[@]}; i++)); do #The length of budget, ${#budget[@]}=2 for budget=(0 20)
    for nballs in ${balls[@]}; do
	cp ${backupdir}/bestrad-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$((i-1))]}).txt bestrad.txt #cp [source] [target]. copy to the current directory

	echo "Running nballs = $nballs, itrial from $((${budget[$((i-1))]}+1)) to ${budget[$i]}, saving to ${backupdir}" #Output information to the command window

	for itrial in $(seq $((${budget[$((i-1))]}+1)) 1 ${budget[$i]}); do
	    if [[ "$problem" = "regpol" && "$InitialGuesses" = "random" ]]; then
    		printf "regpol\n${nballs}\n${itrial}\nrandom\n" | timeout --preserve-status 40 ${PWD}/covering > /dev/null
            elif [[ "$problem" = "regpol" && "$InitialGuesses" = "lattice" ]]; then
		printf "regpol\n${nballs}\n${itrial}\nlattice\nT\nF\n0.03\n0.15\n${nvert}\n" | timeout --preserve-status 40 ${PWD}/covering > /dev/null
	    elif [[ "$problem" = "Agen" && "$InitialGuesses" = "random" ]]; then
	    	#echo "Agen\n${nballs}\n${itrial}\nrandom\n"
		printf "Agen\n${nballs}\n${itrial}\nrandom\n" | timeout --preserve-status 160 ${PWD}/covering > /dev/null
	    elif [[ "$problem" = "Agen" && "$InitialGuesses" = "lattice" ]]; then
	    	#echo "Agen\n${nballs}\n${itrial}\nlattice\nT\nF\n0.03\n0.15\n"
		printf "Agen\n${nballs}\n${itrial}\nlattice\nT\nF\n0.03\n0.15\n" | timeout --preserve-status 160 ${PWD}/covering > /dev/null
		#printf "Agen\n${nballs}\n${itrial}\nT\nF\n0.03\n0.15\n" | timeout --preserve-status 160 ${PWD}/covering > "output_${nballs}_${itrial}.log" 2>&1
	    else
		echo "Usage: $0 Agen|regpol random|lattice [nvert]"s
		exit 1
	    fi

	    # Save iteration information
	    export statdtris=$?
	    if [ -e tabline.txt ] ; then #Check whether a file or directory named tabline.txt exists in the current directory
		sed -e 's/$/  BR/' -i tabline.txt 
		cat tabline.txt >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt #File merging
		bestrad=$(cat bestrad.txt) #Read the content of the file "bestrad.txt" and assign it to the variable "bestrad". 
		bestitrial=$itrial
	    else
		if [ $statdtris -eq 0 ] && [ -e alsolver-interrupted-tabline.txt ]; then 
		    sed -e 's/$/  S/' -i alsolver-interrupted-tabline.txt
		    cat alsolver-interrupted-tabline.txt >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ;
		elif [ $statdtris != 0 ] && [ -e alsolver-interrupted-tabline.txt ]; then 
		    sed -e 's/$/  F/' -i alsolver-interrupted-tabline.txt
		    cat alsolver-interrupted-tabline.txt >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ;
		else
		    printf " no information\n" >> table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ;
		fi
	    fi
	    rm -f tabline.txt alsolver-interrupted-tabline.txt \
	       solver-interrupted-tabline.txt newtkkt-interrupted-tabline.txt
	done

	mv bestrad.txt ${backupdir}/bestrad-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt #Archive the final result to a backup directory
	if [ -e output.dat ]; then
		mv output.dat ${backupdir}/output-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).dat #back up the output.dat after this itrial loop
	fi
	if [ -e picture-solution-final.mp ]; then
		mv picture-solution-final.mp ${backupdir}/picture-solution-final-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).mp #back up the picture-solution-final.mp after this itrial loop
	fi

	if [ "$backupdir" != "." ]; then 
	    if [ "${budget[$i]}" = "0" ]; then 
		mv table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt ${backupdir} #
	    else 
		cat table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt \
		    >> ${backupdir}/table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt 

		rm -f table-${problem}-$(printf "%03d" ${nballs})-$(printf "%06d" ${budget[$i]}).txt 
	    fi
	fi

	# read -p "Press any key to resume ..."
    done
done
