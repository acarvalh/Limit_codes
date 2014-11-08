#!/bin/bash

doBlinding=1

version="v37"
basedir="/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/$version"
limitdirs=("${version}_fitToMgg_nonresSearch_withKinFit" "${version}_fitToMgg_resSearch_withKinFit" "${version}_fitToMgg_resSearch_withRegKinFit" "${version}_fitTo2D_nonresSearch_withKinFit" "${version}_fitTo2D_resSearch_withRegKinFit" "${version}_fitTo2D_resSearch_withKinFit" "${version}_fitToFTR14001_nonresSearch_withKinFit")
doResLimits=("0" "1" "1" "0" "1" "1" "0")
do2DLimits=("0" "0" "0" "1" "1" "1" "1")

#If you only want to run on a subset of directories, edit this array with the appropriate indices.
runLimits=("0" "1" "2" "3" "4" "5" "6")
runLimits=("0")

for i in `echo ${runLimits[@]}`; do

    dir="$basedir/${limitdirs[$i]}"

    if [ ${doResLimits[$i]} == "0" ]; then
	masses=("0")
    else
	masses=("260" "270" "300" "350" "400" "450" "500")
    fi

    if [ ${do2DLimits[$i]} == "0" ]; then
	fitterScript=R2GGBBFitter_mgg_addhiggs.cc
    else
	fitterScript=R2GGBBFitter_2D_addhiggs.cc
    fi

    for imass in `echo ${masses[@]}`; do

	#Prepare the input files
	sed -i "/TString hhiggsggh/c\  TString hhiggsggh = \"$dir/ggh_m125_powheg_8TeV_m${imass}.root\";" $fitterScript
	sed -i "/TString hhiggstth/c\  TString hhiggstth = \"$dir/tth_m125_8TeV_m${imass}.root\";" $fitterScript
	sed -i "/TString hhiggsvbf/c\  TString hhiggsvbf = \"$dir/vbf_m125_8TeV_m${imass}.root\";" $fitterScript
	sed -i "/TString hhiggsvh/c\  TString hhiggsvh = \"$dir/wzh_m125_8TeV_zh_m${imass}.root\";" $fitterScript
	sed -i "/TString hhiggsbbh/c\  TString hhiggsbbh = \"$dir/bbh_m125_8TeV_m${imass}.root\";" $fitterScript
	sed -i "/TString ddata/c\  TString ddata = \"$dir/Data_m${imass}.root\";" $fitterScript

	if [ ${imass} == "0" ]; then
	    sed -i "/TString ssignal/c\  TString ssignal = \"$dir/ggHH_8TeV_m0.root\";" $fitterScript
	elif [ ${imass} == "260" ]; then
	    sed -i "/TString ssignal/c\  TString ssignal = \"$dir/MSSM_m260_8TeV_m260.root\";" $fitterScript
	else
	    sed -i "/TString ssignal/c\  TString ssignal = \"$dir/Radion_m${imass}_8TeV_${imass}.root\";" $fitterScript
	fi

        # the legend
	sed -i -r -e "s/[0-9]+ GeV\"\);/${masses[$j]} GeV\"\);/g" $fitterScript

	echo "Running limits for mass $imass on ${limitdirs[$i]}"

	outputdir="radlim_${limitdirs[$i]}/radlim${imass}"
	mkdir -p $outputdir
	root -l -q runmgg.C >> ${outputdir}/log_radlim${imass}.txt
	mv workspaces/*.root $outputdir
	mv datacards/*.txt $outputdir

       #also colect the plots
	mv databkgoversig* $outputdir
	mv sigmodel* $outputdir
	mv higmodel* $outputdir

        ## create limits root files for each mass
	cd $outputdir
	
	if [ ${doBlinding} == 1 ]    
	then
	    combine -M Asymptotic --run blind hgg.mH125.0_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${imass}_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_higgs.root
	    echo did with Higgs
	        
	    combine -M Asymptotic --run blind hgg.mH125.0_8TeVonecatnohiggs.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${imass}_onecatnohiggs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_onecatnohiggs.root
	    echo did one categ no higgs

	    combine -M Asymptotic --run blind hgg.mH125.0_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.6.mR${imass}_nosyst_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_nosyst_higgs.root
	    echo did with no syst with higgs
	    rm roostats*
	    cd ../..
	else
	    combine -M Asymptotic hgg.mH125.0_8TeV.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${imass}_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_higgs.root
	    echo did with Higgs
	        
	    combine -M Asymptotic hgg.mH125.0_8TeVonecatnohiggs.txt >> higgsCombineTest.Asymptotic.mH125.6.mR${imass}_onecatnohiggs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_onecatnohiggs.root
	    echo did one categ no higgs

	    combine -M Asymptotic hgg.mH125.0_8TeV.txt -S 0 >> higgsCombineTest.Asymptotic.mH125.6.mR${imass}_nosyst_higgs.txt
	    mv higgsCombineTest.Asymptotic.mH120.root higgsCombineTest.Asymptotic.mH125.mR${imass}_nosyst_higgs.root
	    echo did with no syst with higgs
	    rm roostats*
	    cd ../..
	fi
    done # mass
done #dirs


