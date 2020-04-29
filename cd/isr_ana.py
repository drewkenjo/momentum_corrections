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

detname = 'FD' if 'FD' in sys.argv[1] else 'CD'

t0 = datetime.now()

fns = ROOT.vector('string')()
hincs = defaultdict(dict)
for fname in sys.argv[1:]:
  fns.push_back(fname)
  ff = ROOT.TFile(fname)
  for hname in ['hfi','hfiw']:
    for i in range(6):
      if i not in hincs[hname]:
        hincs[hname][i] = ff.Get("{}{}".format(hname,i+1))
      else:
        hincs[hname][i].Add(ff.Get("{}{}".format(hname,i+1)))


df = ROOT.RDataFrame("h22", fns)


hws = []
for i in range(1,7):
  rdf = df.Filter('sec=='+str(i))
  elas = rdf.Filter('thgz>5')
  hws.append(elas.Histo1D(("hw","W at sec "+str(i),100,0,2),"ww"))



def fitpeak(h1,sig,r0,r1):
  minchi2 = float('inf')
  f1 = ROOT.TF1('f1','gaus(0)+pol1(3)',r0,r1)
  f1.SetLineWidth(1)
  mu = 0.9
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
  return pars

musig = []
for i in range(6):
  pars = fitpeak(hws[i], 10*(hws[i].GetBinCenter(2)-hws[i].GetBinCenter(1)), 0.8, 1.05)
  pars = fitpeak(hws[i], pars[2], pars[1]-3*pars[2], pars[1]+3*pars[2])
  musig.extend(pars[1:3])


h1s, h2s = [], []

isrdf = df.Filter('thgz<2')

for i in range(1,7):
  isr = isrdf.Filter('sec=='+str(i))

  h1s.append(isr.Histo1D(('isr/h_W', 'W (sec '+str(i)+');W [GeV]',200,0,5),'ww'))
  h1s.append(isr.Histo1D(('isr/h_phi', 'coplanarity (sec '+')',100,170,190),'fi'))

  h2s.append(isr.Histo2D(('isr/he_thepe', '[electron coverage] ISR: #theta_{e} vs p_{e} (sec '+str(i)+');p_{e} [GeV];#theta_{e}',101,0,11,101,5,40),'pe','the'))
  h2s.append(isr.Histo2D(('isr/hp_thppp', '[proton coverage] ISR: #theta_{p} vs p_{p} (sec '+str(i)+');p_{p} [GeV];#theta_{p}',101,0,3,101,30,75),'pp','thp'))

  h2s.append(isr.Histo2D(('isr/h_e0w', 'ISR: e0 vs W (sec '+str(i)+');W [GeV];E_{0} [GeV]',101,0,5,201,0,12),'ww','e0'))


pdfname = 'isr_hists_'+sys.argv[1].replace('.root','.pdf')

ff = ROOT.TFile('isr_hists_'+sys.argv[1], 'recreate')
ff.cd()

c1 = ROOT.TCanvas("c1",'c1',1100,800)
c1.Print(pdfname+'[')
c1.Divide(3,2,0.0001,0.0001)


for i in range(6):
  c1.cd(i+1)
  hws[i].Draw()
c1.Print(pdfname)


for hname in hincs:
  for i in range(6):
    c1.cd(i+1).SetLogz(1)
    hincs[hname][i].Draw('colz' if 'TH2' in hincs[hname][i].__class__.__name__ else '')
    hincs[hname][i].Write('{}_s{}'.format(hincs[hname][i].GetName(),i+1))
  c1.Print(pdfname)


n1s = len(h1s)/6
for hs in [h1s[i::n1s] for i in range(n1s)]:
  for i in range(6):
    c1.cd(i+1)
    hs[i].Draw()
    hs[i].Write('{}_s{}'.format(hs[i].GetName(),i+1))
  c1.Print(pdfname)

n2s = len(h2s)/6
for hs in [h2s[i::n2s] for i in range(n2s)]:
  for i in range(6):
    c1.cd(i+1).SetLogz(1)

    hs[i].Draw('colz')
    hs[i].Write('{}_s{}'.format(hs[i].GetName(),i+1))
  c1.Print(pdfname)
c1.Print(pdfname+']')

ff.Close()

print((datetime.now()-t0).total_seconds())
