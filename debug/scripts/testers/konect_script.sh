#!/bin/bash

Path="../../input/KONECT-V/KONECT";
log_folder="../../logs/KONECT-OUT-NO-REDUCTION";
mkdir -p $log_folder;

#targets=("1859" "3255" "2300" "1500" "1222" "1345" "1392" "1823" "1668" "1621" "5820" "2438" "5800" "2800" "2747") 

index=0

for i in $Path/* 
do  
    j=$(basename $i);
    date_var=$(date +"%d/%m/%Y")
    time_var=$(date +"%T")

    echo "${j%%.*}," >> ${log_folder}/${j};

    for seed in `seq 1 1`
    do 
        ../../../main -T 60 < $i >> ${log_folder}/${j};
    done
    

    index=$((index+1))
done