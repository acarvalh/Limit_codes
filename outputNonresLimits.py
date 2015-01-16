#!/usr/bin/env python

#run with args <output from combine> <nonres senario> <output file for all limits>
./outputNonresLimits.py $outputdir/$outputfile $limitOutputFile

import sys

if len(sys.argv) < 4:
    echo "Not enough arguments"
    exit 1

    
combineFile=sys.argv[1]
scenario=sys.argv[2]
outputFile=sys.argv[3]



