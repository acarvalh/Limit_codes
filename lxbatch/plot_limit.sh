#!/bin/bash

usage='Usage: -mass <inputMass> -min <inputMin> -max <inpuMax> -step <inputStep> -nsteps <inputNSteps>'



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
      -mass) shift; inputMass=$2;shift;;
      -min) shift; inputMin=$2;shift;;
      -max) shift; inputMax=$2;shift;;
      -step) shift; inputStep=$2;shift;;
      -nsteps) shift; inputNSteps=$2;shift;;
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

if [ "X"${inputStep} == "X" ]
    then
    echo
echo "MISSING INPUT STEP! Please give a valid one!"
    echo $usage
    exit
fi

if [ "X"${inputNSteps} == "X" ]
    then
    echo
echo "MISSING INPUT NUMBER OF STEPS! Please give a valid one!"
    echo $usage
    exit
fi

rm limits.txt

for((i = 0 ; i < ${inputNSteps} ; i++)); do 
 for((j = 0 ; j < ${inputNSteps} ; j++)); do
    mjj_min_tmp=`echo $inputMin-$i*$inputStep|bc| awk '{printf "%f", $0}'`
    mjj_max_tmp=`echo $inputMax+$j*$inputStep|bc| awk '{printf "%f", $0}'`
    mjj_min=$(echo $mjj_min_tmp  | awk ' {sub("\\.*0+$","");print} ')
    mjj_max=$(echo $mjj_max_tmp  | awk ' {sub("\\.*0+$","");print} ')
    echo "Window = ["$mjj_min "," $mjj_max"]"
    for file in `echo "limit_${inputMass}_${mjj_min}_${mjj_max}.txt"`; do
        for exp in `echo "2.5 16.0 50.0 84.0 97.5"`; do
            res=`cat ${file} | 'grep' "Expected" | 'grep' "${exp}%" | cut -d "<" -f 2 `
	    echo "${inputMass} ${mjj_min} ${mjj_max} ${exp} ${res}" >> limits.txt
        done
    done
    for file in `echo "limit_${inputMass}_${mjj_min}_${mjj_max}_cat0.txt"`; do
        for exp in `echo "2.5 16.0 50.0 84.0 97.5"`; do
            res=`cat ${file} | 'grep' "Expected" | 'grep' "${exp}%" | cut -d "<" -f 2 `
	    echo "${inputMass} ${mjj_min} ${mjj_max} ${exp} ${res}" >> limits_cat0.txt
        done
    done
    for file in `echo "limit_${inputMass}_${mjj_min}_${mjj_max}_cat1.txt"`; do
        for exp in `echo "2.5 16.0 50.0 84.0 97.5"`; do
            res=`cat ${file} | 'grep' "Expected" | 'grep' "${exp}%" | cut -d "<" -f 2 `
	    echo "${inputMass} ${mjj_min} ${mjj_max} ${exp} ${res}" >> limits_cat1.txt
        done
    done
 done
done

