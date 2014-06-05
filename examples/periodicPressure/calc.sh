#! /bin/bash

dp=10.0
umax=0.02
outDir=foo

for N in 10 20 40 60 80
do
    for length in 1.0 2.0 4.0
    do
	./periodicPressure $N $dp 0.02 $length $outDir/dp${dp}-N$N-length${length}_ \
	    2> /dev/null \
	    | grep FINAL | cut -d " " -f 2-
    done
done
