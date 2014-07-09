#!/bin/bash

usage='Usage: -r <inputMass> -d <inputMin> -s <inpuMax> -f <inputPoints>'



args=`getopt rd: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi



eval set -- "$args"
for i
  do
  case "$i" in
      -r) shift; inputMass=$2;shift;;
      -d) shift; inputMin=$2;shift;;
      -s) shift; inputMax=$2;shift;;
      -f) shift; inputPoints=$2;shift;;
  esac
done

if [ "X"${inputMass} == "X" ]
    then
    echo
echo "MISSING INPUT MASS! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${inputMin} == "X" ]
    then
    echo
echo "MISSING INPUT MINIMUM! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${inputMax} == "X" ]
    then
    echo
echo "MISSING INPUT MAXIMUM! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${inputPoints} == "X" ]
    then
    echo
echo "MISSING INPUT NUMBER OF POINTS! Please give a valid one!"
    echo $usage
    exit
fi

rm limits.txt

for((i = 0 ; i < ${inputPoints} ; i++)); do 
    mjj_min=$(($inputMin-$i));
    mjj_max=$(($inputMax+$i));
    echo "Mjj window = ["$mjj_min "," $mjj_max"]"
    for file in `echo "limit_${inputMass}_cut_${mjj_min}_${mjj_max}_bkg_reweight.txt"`; do
        for exp in `echo "2.5 16.0 50.0 84.0 97.5"`; do
            res=`cat ${file} | 'grep' "Expected" | 'grep' "${exp}%" | cut -d "<" -f 2 `
	    echo "${inputMass} ${i} ${exp} ${res}" >> limits.txt
        done
    done
done

