#!/bin/bash

Path="../../input/KONECT-V/KONECT";
log_folder="../../logs/KONECT_WITH_MEMORY";
mkdir -p $log_folder;

targets=("1859" "3255" "2300" "1500" "1222" "1345" "1392" "1823" "1668" "1621" "5820" "2438" "5800" "2800" "2747") 

index=0

for i in $Path/* 
do  
    j=$(basename $i);
    date_var=$(date +"%d/%m/%Y")
    time_var=$(date +"%T")

    echo "Testing: $j in $date_var $time_var with target ${targets[$index]}";

    if  [ -e ${log_folder}/${j} ];
    then 
        echo "File does not exist"
    else 
        for seed in `seq 1 1`
        do 
            ../../../main -T 60 < $i >> ${log_folder}/${j};
        done
    fi

    index=$((index+1))
done
    

