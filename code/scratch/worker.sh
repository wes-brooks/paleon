#! /usr/bin/sh

for ((i=1; i<=19 ; i++))
do
    Rscript code/scratch/bootstrap-master.r $i &
done
