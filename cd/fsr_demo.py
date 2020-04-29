#!/usr/bin/python

import sys, ROOT, math
import numpy as np
from array import array
from datetime import datetime
from collections import defaultdict

ROOT.ROOT.EnableImplicitMT(16)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
ROOT.gStyle.SetOptFit(11111)

applycorrection="""
double pars[] = {42.20969552641427, -2.8399875728176522, 0.06307877396588354, -0.0004715997757915601, 102.27538945436321, -6.765941628510161, 0.14921050731600138, -0.001101818195409342, 24.912820827630014, -1.6129034288274415, 0.03460188817226022, -0.00025372562208178013, 56.75698266850059, -3.738750834643781, 0.08227699148433526, -0.0006109981818641908, 161.77315106867073, -10.49066091103833, 0.2269151839670383, -0.0016370803728328394, 299.0233191688738, -19.618606201498572, 0.42813070753183186, -0.003115831496680318};


int ipar = (sec-1)*4;
double thp8 = thp + pars[ipar] + pars[ipar+1]*thp + pars[ipar+2]*thp*thp + pars[ipar+3]*thp*thp*thp;

return thp8;
"""


df = ROOT.RDataFrame("h22", sys.argv[1])\
  .Define('thp8',applycorrection)\
  .Define('e08','0.938*(-1 + 1/(tan(TMath::DegToRad()*the/2)*tan(TMath::DegToRad()*thp8)))')\
  .Filter('thge<3')


h1s = []
for i in range(1,7):
  fsr = df.Filter('sec=='+str(i))
  h1s.append(fsr.Histo1D(("he0","E0 at sec "+str(i),150,8.5,12.5),"e0"))
  h1s.append(fsr.Histo1D(("he08","E0 at sec "+str(i),150,8.5,12.5),"e08"))


def fitpeak0(h1,sig,r0,r1):
  minchi2 = float('inf')
  f1 = ROOT.TF1('f1','gaus(0)+pol1(3)',r0,r1)
  mu = h1.GetBinCenter(h1.GetMaximumBin())
  for sg in np.arange(sig/10,sig,sig/10):
    f1.SetParameters(1,mu,sg,h1.GetMaximum()/10,h1.GetMaximum()/10)
    f1.SetParLimits(1,r0,r1)
    f1.SetParLimits(2,0,r1-r0)
    h1.Fit(f1,"RQ")
    if f1.GetChisquare()<minchi2:
      minchi2 = f1.GetChisquare()
      pars = [f1.GetParameter(i) for i in range(f1.GetNpar())]
  f1.SetParameters(*pars)
  h1.Fit(f1,'RQ')
  return f1


def fitpeak(h1,sig,r0,r1):
  minchi2 = float('inf')
  f1 = ROOT.TF1('f1','[0]*exp(-0.5*((x-[1])/[2])**2) + [0]*[4]*exp(-0.5*((x-[1]-[3])/2/[2])**2) + [0]*[4]*[5]*exp(-0.5*((x-[1]-2*[3])/2/[2])**2)',r0,r1)
  f1.SetLineWidth(1)
  mu = h1.GetBinCenter(h1.GetMaximumBin())
  for sg in np.arange(sig/10,sig,sig/10):
    f1.SetParameters(1,mu,sg,h1.GetMaximum()/10,h1.GetMaximum()/10)
    f1.SetParLimits(1,r0,r1)
    f1.SetParLimits(2,0,r1-r0)
    f1.SetParLimits(3,r0-r1,r1-r0)
    f1.SetParLimits(4,0,1)
    f1.SetParLimits(5,0,1)
    h1.Fit(f1,"RQ")
    if f1.GetChisquare()<minchi2:
      minchi2 = f1.GetChisquare()
      pars = [f1.GetParameter(i) for i in range(f1.GetNpar())]
  f1.SetParameters(*pars)
  h1.Fit(f1,'RQ')
  return f1

conf = {
'he0': {
'def': [9.6,10.8],
5: [9.6, 11.2],
6: [9.4, 11.2]
},
'he08': {
'def': [10.2,11],
5: [10, 11.2],
6: [9.8, 11.5]
},
}

nhs = len(h1s)/6
mgr = defaultdict(lambda: ROOT.TMultiGraph())
for ihist in range(nhs):
  hhs = h1s[ihist::nhs]
  name = 'he0' if 'he0' in hhs[0].GetName() else 'hw'
  gr = ROOT.TGraphErrors()
  gr.SetMarkerSize(2)

  mgr[name].Add(gr)
  igr = mgr[name].GetListOfGraphs().GetSize()

  gr.SetMarkerStyle(20+igr-1)
  gr.SetMarkerColor(igr)
  gr.SetLineColor(igr)

  for i in range(6):
    name = hhs[i].GetName()
    sig = hhs[i].GetBinCenter(10)-hhs[i].GetBinCenter(1)
    f1 = fitpeak(hhs[i], sig, *conf[name].get(i+1, conf[name]['def']))
    mu,sig = f1.GetParameter(1), f1.GetParameter(2)

    gr.SetPoint(gr.GetN(), i+float(ihist)/5, mu)
    gr.SetPointError(gr.GetN()-1, 0, sig)

pdfname = 'fsr_demo_'+sys.argv[1].replace('.root','.pdf')

c1 = ROOT.TCanvas("c1",'c1',1100,800)
c1.Print(pdfname+'[')
c1.Divide(3,2,0.0001,0.0001)

for ihist in range(nhs):
  hhs = h1s[ihist::nhs]
  for i in range(6):
    c1.cd(i+1)
    hhs[i].Draw()
  c1.Print(pdfname)


nom = {'he0': 10.6}
title = {'he0': 'Beam energy via #theta_{e} and #theta_{p};sector number;beam energy [GeV]',
'hw': 'W;sector number;W [GeV]'
}


ll = ROOT.TLine()
ll.SetLineColor(4)
c1.Clear()
for name in mgr:
  leg = ROOT.TLegend(0.8,0.9,1,1)
  mgr[name].Draw("AP")
  mgr[name].SetTitle(title[name])
  leg.AddEntry(mgr[name].GetListOfGraphs()[0],'original momentum','P')
  leg.AddEntry(mgr[name].GetListOfGraphs()[1],'corrected momentum','P')
  leg.Draw()
  ll.DrawLine(mgr[name].GetXaxis().GetXmin(), nom[name], mgr[name].GetXaxis().GetXmax(), nom[name])
  c1.Print(pdfname)
c1.Print(pdfname+']')

