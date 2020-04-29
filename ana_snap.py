#!/apps/bin/python

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

filterstr = "std::array<double, 12> wwms{"+",".join(str(m) for m in musig)+"};return fabs(ww-wwms[sec*2-2])<3*wwms[sec*2-1];"
elasdf = df.Filter(filterstr)
#fsr1 = rdf.Filter('thge<2 && e0>10')
fsr0df = df.Filter('thge>20 || thgz>5')
fsr1df = df.Filter('thge<3')

for i in range(1,7):
  rdf = df.Filter('sec=='+str(i))
  elas = elasdf.Filter('sec=='+str(i))
  fsr0 = fsr0df.Filter('sec=='+str(i))
  fsr1 = fsr1df.Filter('sec=='+str(i))

  h1s.append(rdf.Histo1D(('h_phi{:d}'.format(i),'coplanarity (sec '+')',100,170,190),'fi'))
  h1s.append(rdf.Histo1D(('h_W', 'W (sec '+str(i)+');W [GeV]',200,0,5),'ww'))
  h1s.append(rdf.Histo1D(('h_thgz', '#Theta_{#gamma z} (sec '+str(i)+');#Theta_{#gamma z} [deg]',200,0,10), 'thgz'))
  h1s.append(rdf.Histo1D(('h_thge', '#Theta_{#gamma e} (sec '+str(i)+');#Theta_{#gamma e} [deg]',200,0,40), 'thge'))\

  h2s.append(rdf.Histo2D(('h_e0w', 'e0 vs W (sec '+str(i)+');W [GeV];E_{0} [GeV]',101,0,5,101,0,12),'ww','e0'))
  h2s.append(fsr0.Histo2D(('h_e0wfsr0', 'FSR suppressed: e0 vs W (sec '+str(i)+');W [GeV];E_{0} [GeV]',101,0,5,101,0,12),'ww','e0'))
  h2s.append(fsr1.Histo2D(('fsr/h_e0w', 'FSR: e0 vs W (sec '+str(i)+');W [GeV];E_{0} [GeV]',101,0,5,201,4,15),'ww','e0'))
  h2s.append(fsr1.Histo2D(('fsr/h_thee0', 'FSR: e0 vs #theta_{e} (sec '+str(i)+');#theta_e [#circ];E_{0} [GeV]',101,5,15,101,0,12),'the','e0'))

  h2s.append(rdf.Histo2D(('he_dpe0w', '[electron via E_{0}, #theta_{e}] #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpe0'))
  h2s.append(rdf.Histo2D(('he_dpe1w', '[electron via #theta_{e}, #theta_{p}] #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpe1'))
  h2s.append(fsr0.Histo2D(('he_dpe1wfsr0', '[electron via #theta_{e}, #theta_{p}] FSR suppressed :#Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpe1'))
  h2s.append(fsr1.Histo2D(('fsr/he_dpe1w', '[electron via #theta_{e}, #theta_{p}] FSR: #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpe1'))

  h2s.append(rdf.Histo2D(('hp_dpp0w', '[proton via E_{0}, #theta_{e}] #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpp0'))
  h2s.append(fsr0.Histo2D(('hp_dpp0wfsr0', '[proton via E_{0}, #theta_{e}] FSR suppressed: #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpp0'))
  h2s.append(fsr1.Histo2D(('fsr/hp_dpp0w', '[proton via E_{0}, #theta_{e}] FSR: #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpp0'))
  h2s.append(rdf.Histo2D(('hp_dpp1w', '[proton via #theta_{e}, #theta_{p}] #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpp1'))
  h2s.append(fsr0.Histo2D(('hp_dpp1wfsr0', '[proton via #theta_{e}, #theta_{p}] FSR suppressed: #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpp1'))
  h2s.append(fsr1.Histo2D(('fsr/hp_dpp1w', '[proton via #theta_{e}, #theta_{p}] FSR: #Delta p vs W (sec '+str(i)+');W [GeV];#Delta p [GeV]',101,0,5,101,-1,1),'ww','dpp1'))

  h2s.append(rdf.Histo2D(('h_wpe', 'W vs p (sec '+str(i)+');P_{e} [GeV];W [GeV]',101,0,11, 101,0,5),'pe','ww'))
  h2s.append(rdf.Histo2D(('h_wpp', 'W vs p (sec '+str(i)+');P_{e} [GeV];W [GeV]',101,0,4, 101,0,5),'pp','ww'))

  h2s.append(rdf.Histo2D(('he_dpe0p', '[electron via E_{0}, #theta_{e}] #Delta p vs p_{e} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,8,11,101,-1,1),'pe','dpe0'))
  h2s.append(fsr0.Histo2D(('he_dpe0pfsr0', '[electron via E_{0}, #theta_{e}] FSR suppressed: #Delta p vs p_{e} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,8,11,101,-1,1),'pe','dpe0'))
  h2s.append(fsr1.Histo2D(('fsr/he_dpe0p', '[electron via E_{0}, #theta_{e}] FSR: #Delta p vs p_{e} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,8,11,101,-1,1),'pe','dpe0'))
  h2s.append(fsr1.Histo2D(('fsr/he_dpe0th', '[electron via E_{0}, #theta_{e}] FSR: #Delta p vs #theta_{e} (sec '+str(i)+');#theta_{e} [#circ];#Delta p [GeV]',101,5,15,101,-1,1),'the','dpe0'))

  h2s.append(rdf.Histo2D(('he_dpe1p', '[electron via #theta_{e}, #theta_{p}] #Delta p vs p_{e} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,8,11,101,-1,1),'pe','dpe1'))

  h2s.append(rdf.Histo2D(('hp_dpp0p', '[proton via E_{0}, #theta_{e}] #Delta p vs p_{p} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,0,4,101,-1,1),'pp','dpp0'))
  h2s.append(fsr0.Histo2D(('hp_dpp0pfsr0', '[proton via E_{0}, #theta_{e}] FSR suppressed: #Delta p vs p_{p} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,0,4,101,-1,1),'pp','dpp0'))
  h2s.append(rdf.Histo2D(('hp_dpp1p', '[proton via #theta_{e}, #theta_{p}] #Delta p vs p_{p} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,0,4,101,-1,1),'pp','dpp1'))
  h2s.append(fsr1.Histo2D(('fsr/hp_dpp0pp', '[proton via E_{0}, #theta_{e}] FSR: #Delta p vs p_{p} (sec '+str(i)+');p [GeV];#Delta p [GeV]',101,1,3.5,101,-1,1),'pp','dpp0'))
  h2s.append(fsr1.Histo2D(('fsr/hp_dpp0thp', '[proton via E_{0}, #theta_{e}] FSR: #Delta p vs #theta_{p} (sec '+str(i)+');#theta_{p} [#circ];#Delta p [GeV]',101,30,65,101,-1,1),'thp','dpp0'))

  h2s.append(fsr1.Histo2D(('fsr/hp_dthp0pp', '[proton via E_{0}, #theta_{e}] FSR: #Delta#theta_{p} vs p_{p} (sec '+str(i)+');p_{p} [GeV];#Delta#theta_{p}',101,1,3.5,101,-10,10),'pp','dthp0'))

  h2s.append(rdf.Histo2D(('hp_dthp0thp', '[proton via E_{0}, #theta_{e}] #Delta#theta_{p} vs #theta_{p} (sec '+str(i)+');#theta_{p};#Delta#theta_{p}',101,30,65,101,-10,10),'thp','dthp0'))
  h2s.append(fsr1.Histo2D(('fsr/hp_dthp0thp', '[proton via E_{0}, #theta_{e}] FSR: #Delta#theta_{p} vs #theta_{p} (sec '+str(i)+');#theta_{p};#Delta#theta_{p}',101,30,65,101,-10,10),'thp','dthp0'))

  #h2s.append(rdf.Histo2D(('h_pe0pe1', 'p vs p (sec '+str(i)+');p [GeV];p [GeV]',101,0,11,101,0,11),'pe','pe1'))
  #h2s.append(rdf.Histo2D(('h_thgethgz', ';e#gamma angle;Z#gamma angle',101,0,40,101,0,40),'thge','thgz'))

  h2s.append(fsr1.Histo2D(('fsr/hp_thppp', '[proton coverage] FSR: #theta_{p} vs p_{p} (sec '+str(i)+');p_{p} [GeV];#theta_{p}',101,1,3.5,101,30,65),'pp','thp'))

  h2s.append(rdf.Histo2D(('h_thgew', '#Theta_{#gamma e} vs W (sec '+str(i)+');W [GeV];#Theta_{#gamma e} [deg]',200,0,5,200,0,40), 'ww','thge'))
  h2s.append(rdf.Histo2D(('h_thgzw', '#Theta_{#gamma z} vs W (sec '+str(i)+');W [GeV];#Theta_{#gamma z} [deg]',200,0,5,200,0,40), 'ww','thgz'))

  h1s.append(elas.Histo1D(('elas/h_W', 'W (sec '+str(i)+');W [GeV]',200,0,2),'ww'))
  h2s.append(elas.Histo2D(('elas/he_thepe', '[electron] #theta_{e} vs p_{e} (sec '+str(i)+');p_{e} [GeV];#theta_{e} [#circ]',101,8,10.5,101,5,15),'pe','the'))
  h2s.append(elas.Histo2D(('elas/he_dpe0pe', '[electron via E_{0}, #theta_{e}] #Delta p_{e} vs p_{e} (sec '+str(i)+');p_{e} [GeV];#Delta p_{e} [GeV]',101,8,10.5,101,-0.25,0.25),'pe','dpe0'))
  h2s.append(elas.Histo2D(('elas/he_dpe0the', '[electron via E_{0}, #theta_{e}] #Delta p_{e} vs #theta_{e} (sec '+str(i)+');#theta_{e} [GeV];#Delta p_{e} [GeV]',101,5,14.5,101,-0.25,0.25),'the','dpe0'))


name = 'all_'+sys.argv[1].replace('.root','.pdf')

ff = ROOT.TFile(sys.argv[1].replace('snap','hists'), 'recreate')
ff.cd()

c1 = ROOT.TCanvas("c1",'c1',1100,800)
c1.Print(name+'[')
c1.Divide(3,2,0.0001,0.0001)

hnames = list(set([h1.GetName() for h1 in hws + h1s + h2s]))
hdict = {hname:[h1 for h1 in hws + h1s + h2s if h1.GetName()==hname]  for hname in hnames}
hdict.update(hincs)



for i in range(6):
  c1.cd(i+1)
  hws[i].Draw()
  f1 = hws[i].GetFunction('f1')
c1.Print(name)


for hname in hincs:
  for i in range(6):
    c1.cd(i+1).SetLogz(1)
    hincs[hname][i].Draw('colz' if 'TH2' in hincs[hname][i].__class__.__name__ else '')
    hincs[hname][i].Write('{}_s{}'.format(hincs[hname][i].GetName(),i))
  c1.Print(name)

n1s = len(h1s)/6
for hs in [h1s[i::n1s] for i in range(n1s)]:
  for i in range(6):
    c1.cd(i+1)
    hs[i].Draw()
    hs[i].Write('{}_s{}'.format(hs[i].GetName(),i))
  c1.Print(name)

ll = ROOT.TLine()
n2s = len(h2s)/6
for hs in [h2s[i::n2s] for i in range(n2s)]:
  c1.Clear()
  c1.Divide(3,2, 0.00001,0.00001)

  for i in range(6):
    c1.cd(i+1).SetLogz(1)
    ROOT.gPad.SetLeftMargin(0.15)

    hs[i].Draw('colz')
    hs[i].Write('{}_s{}'.format(hs[i].GetName(),i))
    if 'he0w' in hs[i].GetName():
      ll.DrawLine(hs[i].GetXaxis().GetXmin(),10.6,hs[i].GetXaxis().GetXmax(),10.6)

  c1.Print(name)

c1.Print(name+']')

ff.Close()

print((datetime.now()-t0).total_seconds())
