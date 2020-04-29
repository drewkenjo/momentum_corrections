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
double pars[] = {0.9659327216385852, -0.31326066844446493, 0.03354932668101431, -0.001211903430960084, 0.18624258921255368, -0.04582791746250911, 0.0036262968892478604, -0.00011017970070603751, -0.7906105326209454, 0.2633876853743531, -0.031063980087627507, 0.0011977810376170585, 0.03382326048166492, -0.056685006505159126, 0.010093793724685282, -0.0004776750705265613, -0.47102396528246576, 0.12100014919041127, -0.011120646412637672, 0.00036102119874576057, -1.3001077988561707, 0.3951973135970893, -0.040136029638208295, 0.0013359750635126284};

int ipar = (sec-1)*4;
double pe9 = pe + pars[ipar] + pars[ipar+1]*the + pars[ipar+2]*the*the + pars[ipar+3]*the*the*the;

//if(rdfentry_<100) std::cout<<the<<" "<<pe9<<std::endl;

return pe9;
"""

E0 = 10.6
if '7GeV' in sys.argv[1]:
  E0 = 7.54626
elif '6GeV' in sys.argv[1]:
  E0 = 6.53536


wcalc="""
double E0=10.6, Mpro=0.938;

TLorentzVector beam(0,0,E0,E0);
TLorentzVector targ(0,0,0,Mpro);
TLorentzVector ele;
TVector3 evec(ex,ey,ez);
evec.SetMag(pe9);
ele.SetVectM(evec, 0);

TLorentzVector wvec = beam+targ-ele;
return wvec.M();
"""


df = ROOT.RDataFrame("h22", sys.argv[1])\
  .Define('pe9',applycorrection)\
  .Define('ww9',wcalc)\
  .Define('e02',"(0.938*pe)/(0.938 - pe + pe*cos(the*TMath::DegToRad()))")\
  .Define('e09',"(0.938*pe9)/(0.938 - pe9 + pe9*cos(the*TMath::DegToRad()))")\


h1s = []
for i in range(1,7):
  rdf = df.Filter('sec=='+str(i))
  elas = rdf.Filter('thgz>5')
  h1s.append(elas.Histo1D(("hw","W at sec "+str(i),100,0,2),"ww"))
  h1s.append(elas.Histo1D(("hw9","W at sec "+str(i),100,0,2),"ww9"))
  h1s.append(elas.Histo1D(("he0","E0 at sec "+str(i),100,9.5,11.5),"e0"))
  h1s.append(elas.Histo1D(("he02","E0 at sec "+str(i),100,9.5,11.5),"e02"))
  h1s.append(elas.Histo1D(("he09","E0 at sec "+str(i),100,9.5,11.5),"e09"))


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
'def': [10,10.8],
5: [10,11],
6: [9.8,10.8],
},
'he02': {
'def': [10.4,10.8]
},
'he09': {
'def': [10.4,10.8]
},
'hw': {
'def': [0.7,1.1]
},
'hw9': {
'def': [0.7,1.1]
}
}

nhs = len(h1s)/6
grs = {}
for ihist in range(nhs):
  hhs = h1s[ihist::nhs]
  name = hhs[0].GetName()
  gr = ROOT.TGraphErrors()
  gr.SetMarkerSize(2)
  grs[hhs[0].GetName()] = gr
  for i in range(6):
    name = hhs[i].GetName()
    sig = hhs[i].GetBinCenter(10)-hhs[i].GetBinCenter(1)
    f1 = fitpeak(hhs[i], sig, *conf[name].get(i+1, conf[name]['def']))
    mu,sig = f1.GetParameter(1), f1.GetParameter(2)

    gr.SetPoint(gr.GetN(), i+float(ihist)/5, mu)
    gr.SetPointError(gr.GetN()-1, 0, sig)


pdfname = 'elas_demo_'+sys.argv[1].replace('.root','.pdf')

c1 = ROOT.TCanvas("c1",'c1',1100,800)
c1.Print(pdfname+'[')
c1.Divide(3,2,0.0001,0.0001)

ll = ROOT.TLine()
for ihist in range(nhs):
  hhs = h1s[ihist::nhs]
#  f1s = []
  for i in range(6):
    c1.cd(i+1)
    hhs[i].Draw()

    ymax = hhs[i].GetMaximum()*1.05
    xline = 0.938 if 'hw' in hhs[i].GetName() else 10.6
    ll.DrawLine(xline,0,xline,ymax)

#    f1 = hhs[i].GetListOfFunctions()[0]
#    pars = [f1.GetParameter(i) for i in range(f1.GetNpar())]

#    f2 = ROOT.TF1("f2",'gaus(0)+gaus(3)',f1.GetXaxis().GetXmin(),f1.GetXaxis().GetXmax())
#    f2.SetParameters(pars[0]*pars[4], pars[1]+pars[3], pars[2]*2, pars[0]*pars[4]*pars[5], pars[1]+2*pars[3], pars[2]*2)
#    f2.SetLineColor(1)
#    f2.Draw("same")

#    f1s.append(f2)
  c1.Print(pdfname)



nom = {'he0': 10.6, 'he02': 10.6, 'hw': 0.938}
title = {'he02': 'Beam energy via P_{e} and #theta_{e};sector number;beam energy [GeV]',
'he0': 'Beam energy;sector number;beam energy [GeV]',
'hw': 'W;sector number;W [GeV]'
}
legend = {'he02': 'via P_{e} and #theta_{e}',
'he09': 'via corrected P_{e} and #theta_{e}',
'he0': 'via #theta_{e} and #theta_{p}',
'hw': 'original W',
'hw9': 'corrected W',
}


ll = ROOT.TLine()
ll.SetLineColor(4)
c1.Clear()
for names in [['he0','he02'], ['he02','he09'], ['hw','hw9']]:
  leg = ROOT.TLegend(0.8,0.9,1,1)
  mgr = ROOT.TMultiGraph()
  for nm,igr in zip(names,range(len(names))):
    grs[nm].SetMarkerStyle(20+igr)
    grs[nm].SetMarkerColor(igr+1)
    grs[nm].SetLineColor(igr+1)
    leg.AddEntry(grs[nm], legend[nm], 'P')
    mgr.Add(grs[nm].Clone())
  mgr.Draw("AP")
  mgr.SetTitle(title[names[0]])
  leg.Draw()
  ll.DrawLine(mgr.GetXaxis().GetXmin(), nom[names[0]], mgr.GetXaxis().GetXmax(), nom[names[0]])
  c1.Print(pdfname)

c1.Print(pdfname+']')

