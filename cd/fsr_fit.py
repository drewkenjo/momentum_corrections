#!/usr/bin/python

import sys, math, ROOT
import numpy as np
from collections import defaultdict
from argparse import Namespace

ROOT.TH1.AddDirectory(0)

keys = defaultdict(list)
ff = ROOT.TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
  if 'fsr/' in kk.GetName():
    keys[kk.ReadObj().GetName()].append((kk.GetName(), kk.ReadObj()))
for kk in keys:
  keys[kk].sort()

hhs = { kk:[tpl[1] for tpl in keys[kk]] for kk in keys}

#n = Namespace(**hhs)

fname = 'fsr_fit_'+sys.argv[1].replace('.root','.pdf')

c1 = ROOT.TCanvas("c1","c1",1100,800)
c1.Divide(3,2,0.0001,0.0001)
c1.Print(fname+'[')
for name in sorted(hhs):
  for i in range(6):
    c1.cd(i+1).SetLogz()
    hhs[name][i].Draw("colz")
  c1.Print(fname)




def fitpeak(h1,sig,r0,r1):
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




fitlims = [
{'def':[-2,2], },
{'def':[-2,2], },
{'def':[-2,2], },
{'def':[-2,2], },
{'def':[-2,2], },
{'def':[-3.2,2], },
]


xlims = np.arange(40,55,1.5)

def fitslices(h2s):
  slices = []
  for h2,isec in zip(h2s,range(len(h2s))):
    gr = ROOT.TGraphErrors()
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1)
    hx = h2.ProjectionX("hx")
    xbins = [hx.GetXaxis().FindBin(i) for i in xlims]
    slices.append({'fits':[]})
    for ix0,ix1 in zip(xbins[:-1],xbins[1:]):
      hy = h2.ProjectionY("hdp{}".format(ix0),ix0,ix1)
      f1 = fitpeak(hy, hy.GetBinCenter(10)-hy.GetBinCenter(1), *fitlims[isec].get(ix0,fitlims[isec]['def']))
      slices[-1]['fits'].append(hy)
      hx.GetXaxis().SetRange(ix0,ix1)
      gr.SetPoint(gr.GetN(), hx.GetMean(), f1.GetParameter(1))
      gr.SetPointError(gr.GetN()-1, 0, f1.GetParError(1))
    slices[-1]['gr'] = gr
    slices[-1]['h2'] = h2
  return slices



for i in range(6):
  c1.cd(i+1).SetLogz(0)
  hhs['fsr/hp_dthp0thp'][i].Draw("colz")
c1.Print(fname)


parlist = []
slices = fitslices(hhs['fsr/hp_dthp0thp'])

for i in range(6):
  c1.cd(i+1).SetLogz()
  slices[i]['gr'].Fit('pol3','Q')
  slices[i]['h2'].Draw('colz')
  slices[i]['gr'].Draw("P")
  f1 = slices[i]['gr'].GetListOfFunctions()[0]
  parlist.extend([f1.GetParameter(i) for i in range(f1.GetNpar())])
c1.Print(fname)

print(parlist)

for secs in slices:
  fits = secs['fits']
  hnum = int(round(math.sqrt(len(fits))))
  vnum = int(math.ceil(float(len(fits))/hnum))
  c1.Clear()
  c1.Divide(hnum,vnum)
  ipad = 0
  for hb in fits:
    ipad += 1
    c1.cd(ipad)
    hb.Draw()
  c1.Print(fname)
  

c1.Print(fname+']')

ff.Close()
