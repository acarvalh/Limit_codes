#!/usr/bin/env python
import string, sys, os, getopt, subprocess, time, shutil
from math import log10, pow
import numpy
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1F, TH2F, TLatex, TLegend, TStyle, gStyle, TMultiGraph, TGraph, TGraphErrors
import math as m

def usage():
    print "Usage: plot_limit_1edge.py --inputfile=[inputfile] --mass=[mass] --min=[mjjmin] --max=[mjjmax] --step=[step] --nsteps=[nsteps]"
    
try:
     opts, args = getopt.getopt(sys.argv[1:], "m:d:r:g:so:", ["inputfile=","mass=","min=","max=","step=","nsteps="])

except getopt.GetoptError:
     #* print help information and exit:*
     usage()
     sys.exit(2)

mass = ''
mjjmin = ''
mjjmax = ''

for opt, arg in opts:
    
     if opt in ("--inputfile"):
        inputfile = arg
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
 
file = inputfile
name  = "cut_Mggjj_study_reweight_m" + mass # + "_zoom"
title = "Expected Limits (normalized to the best value) at mass = " + mass + " GeV in both categories. Use of all bkg samples"
XaxisTittle = "Mggjj_Min (GeV)"
YaxisTittle = "Expected Limit (normalized to the best value)"

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

h1 = TH1F("h1",title,int(nsteps),float(mjjmin)-float(step)*(int(nsteps)-1)-float(step)/2,float(mjjmin)+float(step)/2)
h1.SetStats(0)

x = numpy.asarray(x, dtype='float')
l_50 = numpy.asarray(l_50, dtype='float')
l_50_sort = numpy.asarray(l_50_sort, dtype='float')
l_50_sort.sort()

norm=l_50_sort[0]
l_50 = [i/norm for i in l_50]

for i in xrange(int(nsteps)) :
   #print i," - ",x[i]," - ",y[i]," - ",l_50[i]
   #if l_50[i] > 3.: 
   #   continue   
   h1.SetBinContent(h1.FindBin(x[i]),l_50[i])


h1.GetXaxis().SetTitle(XaxisTittle)
h1.GetYaxis().SetTitle(YaxisTittle)
h1.GetYaxis().SetRangeUser(0.99,1.2)
h1.SetMarkerStyle(20)
h1.SetMarkerSize(0.5)
h1.Draw("P")

c1.Update()

c1.Print(name + ".pdf");
c1.Print(name + ".gif");
c1.Print(name + ".root");


