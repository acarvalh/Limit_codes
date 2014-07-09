#!/usr/bin/env python
import string, sys, os, getopt, subprocess, time, shutil
from math import log10, pow
import numpy
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TLatex, TLegend, TStyle, gStyle, TMultiGraph, TGraph, TGraphErrors
import math as m

def usage():
    print "Usage: plot_limit.py --mass=[mass] --min=[mjjmin] --max=[mjjmax]"
    
try:
     opts, args = getopt.getopt(sys.argv[1:], "m:d:r:g:so:", ["mass=","min=","max="])

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
 
file = "limits.txt"
name  = "cut_Mjj_study_reweight_m" + mass # + "_zoom"
title = "Expected Limits at mass = " + mass + " GeV in both category. Use of all bkg samples"
XaxisTittle = "Mjj = [" + mjjmin + "-i," + mjjmax + "+i] GeV"
YaxisTittle = "Expected Limit"

c1 = TCanvas()
c1.SetGridy()

x = []
l_2_5 =[]
l_16 =[]
l_50 =[]
l_84 =[]
l_97_5 =[]

with open(file) as data:
    for line in data:
        mass_, cut1_, exp_, res_ = line.split()
        if str(2.5) in exp_ :
            l_2_5.append(float(res_))
            x.append(float(cut1_))
        elif str(16.0) in exp_:
            l_16.append(float(res_))
        elif str(50.0) in exp_:
            l_50.append(float(res_))
        elif str(84.0) in exp_:
            l_84.append(float(res_))
        elif str(97.5) in exp_:
            l_97_5.append(float(res_))

norm = 1.
g = TMultiGraph()
g.SetTitle(title)
g2 = TMultiGraph()


x = numpy.asarray(x, dtype='float')
x_e = []
for i in range(len(x)):
    x_e.append(float(0.02))
x_e = numpy.asarray(x_e, dtype='float')


l_2_5 = [i/norm for i in l_2_5]
l_16 = [i/norm for i in l_16]
l_50 = [i/norm for i in l_50]
l_84 = [i/norm for i in l_84]
l_97_5 = [i/norm for i in l_97_5]

l_2_5 = numpy.asarray(l_2_5, dtype='float')
l_2_5_m = numpy.asarray((l_2_5 + l_16)/2, dtype='float')
l_2_5_e = numpy.asarray(-(l_2_5 - l_16)/2, dtype='float')
g_l_2_5= TGraphErrors(len(l_2_5), x, l_2_5_m, x_e, l_2_5_e)
g_l_2_5.SetFillColor(5)
g_l_2_5.SetFillStyle(1001)
g.Add(g_l_2_5)

l_16 = numpy.asarray(l_16, dtype='float')
l_16_m = numpy.asarray((l_16+l_50)/2, dtype='float')
l_16_e = numpy.asarray(-(l_16-l_50)/2, dtype='float')
g_l_16= TGraphErrors(len(l_16), x, l_16_m, x_e, l_16_e)
g_l_16.SetFillColor(3)
g_l_16.SetFillStyle(1001)
g.Add(g_l_16)

l_50 = numpy.asarray(l_50, dtype='float')
g_l_50= TGraphErrors(len(l_50), x, l_50)
g_l_50.SetMarkerColor(1)
g_l_50.SetMarkerStyle(2)
g_l_50.SetMarkerSize(1.5)
g2.Add(g_l_50)


l_84 = numpy.asarray(l_84, dtype='float')
l_84_m = numpy.asarray((l_84+l_50)/2, dtype='float')
l_84_e = numpy.asarray(-(l_84-l_50)/2, dtype='float')
g_l_84= TGraphErrors(len(l_84), x, l_84_m, x_e, l_84_e)
g_l_84.SetFillColor(3)
g_l_84.SetFillStyle(1001)
g.Add(g_l_84)

l_97_5 = numpy.asarray(l_97_5, dtype='float')
l_97_5_m = numpy.asarray((l_97_5+l_84)/2, dtype='float')
l_97_5_e = numpy.asarray(-(l_97_5-l_84)/2, dtype='float')
g_l_97_5= TGraphErrors(len(l_97_5), x, l_97_5_m, x_e, l_97_5_e)
g_l_97_5.SetFillColor(5)
g_l_97_5.SetFillStyle(1001)
g.Add(g_l_97_5)
#g.SetMaximum(1.05 * max(l_97_5))
#g.SetMaximum(2.5)
#g.SetMaximum(1.05)
#g.SetMinimum(.4)
#g.SetMinimum(.85)
g.Draw("a2")
g2.Draw("f")
g.GetXaxis().SetTitle(XaxisTittle)
g.GetYaxis().SetTitle(YaxisTittle)


c1.Update()

c1.Print(name + ".pdf");
c1.Print(name + ".gif");
c1.Print(name + ".root");


