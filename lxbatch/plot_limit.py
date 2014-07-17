#!/usr/bin/env python
import string, sys, os, getopt, subprocess, time, shutil
from math import log10, pow
import numpy
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TH2F, TLatex, TLegend, TStyle, gStyle, TMultiGraph, TGraph, TGraphErrors
import math as m

def usage():
    print "Usage: plot_limit.py --mass=[mass] --min=[mjjmin] --max=[mjjmax] --step=[step] --nsteps=[nsteps]"
    
try:
     opts, args = getopt.getopt(sys.argv[1:], "m:d:r:g:so:", ["mass=","min=","max=","step=","nsteps="])

except getopt.GetoptError:
     #* print help information and exit:*
     usage()
     sys.exit(2)

mass = ''
mjjmin = ''
mjjmax = ''

for opt, arg in opts:
    
     if opt in ("--mass"):
        mass = arg
     if opt in ("--min"):
        mjjmin= arg     
     if opt in ("--max"):
        mjjmax= arg       
     if opt in ("--step"):
        step= arg    
     if opt in ("--nsteps"):
        nsteps= arg    

print "Mass    = ",mass
print "Min     = ",mjjmin 
print "Max     = ",mjjmax
print "Step    = ",step
print "nSteps  = ",nsteps
 
file = "limits.txt"
name  = "cut_Mjj_study_reweight_m" + mass # + "_zoom"
title = "Expected Limits at mass = " + mass + " GeV in both category. Use of all bkg samples"
XaxisTittle = "Mjj_Min (GeV)"
YaxisTittle = "Mjj_Max (GeV)"
#ZaxisTittle = "Expected Limit"

c1 = TCanvas()

x = []
y = []
l_50 =[]
l_50_sort =[]

with open(file) as data:
    for line in data:
        mass_, cut1_, cut2_, exp_, res_ = line.split()
        if str(50.0) in exp_ :
            l_50.append(float(res_))
            l_50_sort.append(float(res_))
            x.append(float(cut1_))
            y.append(float(cut2_))

h2 = TH2F("h2",title,int(nsteps),float(mjjmin)-float(step)*(int(nsteps)-1)-float(step)/2,float(mjjmin)+float(step)/2,int(nsteps),float(mjjmax)-float(step)/2,float(mjjmax)+float(step)*(int(nsteps)-1)+float(step)/2)
h2.SetStats(0)

x = numpy.asarray(x, dtype='float')
y = numpy.asarray(y, dtype='float')
l_50 = numpy.asarray(l_50, dtype='float')
l_50_sort = numpy.asarray(l_50_sort, dtype='float')
l_50_sort.sort()

norm=l_50_sort[0]
l_50 = [i/norm for i in l_50]

for i in xrange(int(nsteps)*int(nsteps)) :
   #print i," - ",x[i]," - ",y[i]," - ",l_50[i]
   #if l_50[i] > 3.: 
   #   continue   
   h2.SetBinContent(h2.FindBin(x[i],y[i]),l_50[i])


h2.GetXaxis().SetTitle(XaxisTittle)
h2.GetYaxis().SetTitle(YaxisTittle)
#h2.GetZaxis().SetTitle(ZaxisTittle)
h2.GetZaxis().SetRangeUser(0.95,1.05)
h2.Draw("colz")

c1.Update()

c1.Print(name + ".pdf");
c1.Print(name + ".gif");
c1.Print(name + ".root");


