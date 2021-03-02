#!/bin/bash

while getopts a:r:t:i: option
	do
		case "${option}"
		in
		a) assembly=${OPTARG};;
		r) reads=${OPTARG};;
		t) threads=${OPTARG};;
		i) iteration_num=${OPTARG};;
		esac
	done

#echo ${assembly} ${reads} ${threads} ${iteration_num}

# assembly=
# reads=
# threads=4
# iteration_num=4

iteration=1
while [ $iteration -le $iteration_num ]

	do
	echo Iteration: ${iteration}
	minimap2 -ax map-ont -t ${threads} ${assembly} ${reads} | gzip -1 > minimap2.${iteration}.sam.gz
	/home/fay-wei/bin/racon/build/bin/racon -m 8 -x -6 -g -8 -w 500 -t ${threads} ${reads} minimap2.${iteration}.sam.gz ${assembly} > racon.${iteration}.fa

	assembly=racon.${iteration}.fa
	iteration=$((iteration+1))

	done
