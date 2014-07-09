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

sed -i -r -e "s/runfits\([0-9]+/runfits\(${inputMass}/g" runfits.C 
sed -i -r -e "s/m[0-9]+.root/m${inputMass}.root/g" R2GGBBFitter_mtot_range.cc 
sed -i -r -e "s/Radion_m[0-9]+_8TeV_m[0-9]+.root/Radion_m${inputMass}_8TeV_m${inputMass}.root/g" R2GGBBFitter_mgg_addhiggs.cc

#windows for cuts
for((i = 0 ; i < ${inputPoints} ; i++)); do 
      mjj_min=$(($inputMin-$i))
      mjj_max=$(($inputMax+$i))
      echo "Mjj window = ["$mjj_min "," $mjj_max"]"	
      sed -i -r -e "s/mjj > [0-9]+ && mjj < [0-9]+/mjj > ${mjj_min} \&\& mjj < ${mjj_max}/g" R2GGBBFitter_mtot_range.cc
      root -l -q -b runfits.C
      cp workspaces/*.root .
      mv datacards/hgg.mH${inputMass}.0_8TeVrep.txt datacards/hgg.mH${inputMass}_${mjj_min}_${mjj_max}.0_8TeVrep.txt
      combine -M Asymptotic datacards/hgg.mH${inputMass}_${mjj_min}_${mjj_max}.0_8TeVrep.txt > limit_${inputMass}_cut_${mjj_min}_${mjj_max}_bkg_reweight.txt
done
