#!/bin/bash

if [ $# -eq 0 ]
then
    echo "Usage: No files from chipatlas supplied"
    exit 1
fi 

for i in $*
do
    echo $i
    name=${i:0:14}
    tail=${i:24:33}
    echo $name$tail
    sed '1d' $i |\
    awk -F '[\t;=%]' '{print $1"\t"$2"\t"$3"\t"$7}'  > $name$tail 
    ## Currently hardcoded to output the files to data folder the filtered files 
done