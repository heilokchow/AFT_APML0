#!/bin/bash

for c in 0
do
	for n in 100 200 300
	do
		for p in 500 2000 10000
		do
			make simulation N_SIM=$n P_SIM=$p N_REP=100 N_THREAD=1 CV_STRUCT=$c
			mkdir $c-$n-$p
			mv test1 ./$c-$n-$p
			cd $c-$n-$p
			./test1 &
			cd ..
		done
	done
done

for c in 1 2
do
	for n in 200 400
	do
		for p in 500 2000 10000
		do
			make simulation N_SIM=$n P_SIM=$p N_REP=100 N_THREAD=1 CV_STRUCT=$c
			mkdir $c-$n-$p
			mv test1 ./$c-$n-$p
			cd $c-$n-$p
			./test1 &
			cd ..
		done
	done
done
