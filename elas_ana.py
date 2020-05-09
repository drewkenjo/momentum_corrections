#!/usr/bin/python

import sys, ROOT, math
import numpy as np
from array import array
from datetime import datetime
from collections import defaultdict
from Conf import HistConf

ROOT.ROOT.EnableImplicitMT(16)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
ROOT.gStyle.SetOptFit(11111)

t0 = datetime.now()

fns = ROOT.vector('string')()
hincs = defaultdict(dict)
for fname in sys.argv[1:]:
  fns.push_back(fname)
  ff = ROOT.TFile(fname)
  E0 = str(ff.Get("E0"))
  for hname in ['hfi','hfiw']:
    for i in range(6):
      if i not in hincs[hname]:
        hincs[hname][i] = ff.Get("{}_s{}".format(hname,i+1))
      else:
        hincs[hname][i].Add(ff.Get("{}_s{}".format(hname,i+1)))


df = ROOT.RDataFrame("h22", fns)


hws = []
for i in range(1,7):
  rdf = df.Filter('sec=='+str(i))
  elas = rdf.Filter('thgz>5')
  hws.append(elas.Histo1D(("elas/h_Wfit","W fit",100,0,2),"ww"))



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


hhs = []

fname = sys.argv[1].split('/')[-1]
hconf= HistConf(fname)
hlist = [cn for cn in hconf.histos if cn[0][0].startswith('elas/')]


thp8correction="double pars[] = {" + ",".join(str(vv) for vv in hconf.thp8pars) + """};
int ipar = (sec-1)*4;
double thp8 = thp + pars[ipar] + pars[ipar+1]*thp + pars[ipar+2]*thp*thp + pars[ipar+3]*thp*thp*thp;

return thp8;
"""

pe9correction="double pars[] = {" + ",".join(str(vv) for vv in hconf.pe9pars) + """};
int ipar = (sec-1)*4;
double pe9 = pe + pars[ipar] + pars[ipar+1]*the + pars[ipar+2]*the*the + pars[ipar+3]*the*the*the;

//if(rdfentry_<100) std::cout<<the<<" "<<pe9<<std::endl;

return pe9;
"""

wcalc = "double E0=" + E0 + ", Mpro=0.938;" + """
TLorentzVector beam(0,0,E0,E0);
TLorentzVector targ(0,0,0,Mpro);
TLorentzVector ele;
TVector3 evec(ex,ey,ez);
evec.SetMag(pe9);
ele.SetVectM(evec, 0);

TLorentzVector wvec = beam+targ-ele;
return wvec.M();
"""



filterstr = "std::array<double, 12> wwms{"+",".join(str(m) for m in musig)+"};return fabs(ww-wwms[sec*2-2])<3*wwms[sec*2-1];"
elasdf = df.Filter(filterstr)\
  .Define('e02', "(0.938*pe)/(0.938 - pe + pe*cos(the*TMath::DegToRad()))")\
  .Define('e03', "double cthe = 1+0.938*(1/"+E0+" - 1/pe); return 0.938*(-1+1/(sqrt((1-cthe)/(1+cthe))*tan(thp*TMath::DegToRad())));")\
  .Define('thp8', thp8correction)\
  .Define('pe9', pe9correction)\
  .Define('ww9', wcalc)\
  .Define('e09', "(0.938*pe9)/(0.938 - pe9 + pe9*cos(the*TMath::DegToRad()))")\
  .Define('e08', "0.938*(-1+1/(tan(the*TMath::DegToRad()/2)*tan(thp8*TMath::DegToRad())))")\



for i in range(1,7):
  elas =  elasdf.Filter('sec=='+str(i))

  for cnf in hlist:
    if len(cnf)==2:
      hhs.append(elas.Histo1D(*cnf))
    else:
      hhs.append(elas.Histo2D(*cnf))

pdfname = 'elas_hists_'+fname.replace('.root','.pdf')

#ff = ROOT.TFile('elas_hists_'+sys.argv[1], 'recreate')
ff = ROOT.TFile(sys.argv[1], 'update')
ff.cd()

c1 = ROOT.TCanvas("c1",'c1',1100,800)
c1.Print(pdfname+'[')
c1.Divide(3,2,0.0001,0.0001)


ROOT.TObjString(filterstr).Write("elasticfilter", ROOT.TObject.kOverwrite)

for hname in hincs:
  for i in range(6):
    c1.cd(i+1).SetLogz(1)
    hincs[hname][i].Draw('colz' if 'TH2' in hincs[hname][i].__class__.__name__ else '')
    #hincs[hname][i].Write('{}_s{}'.format(hincs[hname][i].GetName(),i+1))
  c1.Print(pdfname)


n1s = len(hhs)/6
for hs in [hws] + [hhs[i::n1s] for i in range(n1s)]:
  for i in range(6):
    is2d = 'TH2' in hs[i].__class__.__name__
    c1.cd(i+1).SetLogz(1)
    hs[i].SetTitle(hs[i].GetTitle()+" (sec "+str(i+1)+")")
    hs[i].Draw('colz' if is2d else '')
    hs[i].Write('{}_s{}'.format(hs[i].GetName(),i+1), ROOT.TObject.kOverwrite)
  c1.Print(pdfname)

c1.Print(pdfname+']')

ff.Close()

print((datetime.now()-t0).total_seconds())
