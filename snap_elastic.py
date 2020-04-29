#!/apps/bin/python

import sys, ROOT, math, re
from array import array

ROOT.ROOT.EnableImplicitMT(12)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(0)
ROOT.gStyle.SetOptFit(11111)

idet,detname = [2,'FD'] if 'FD' in sys.argv else [4,'CD']
fnames = [fn for fn in sys.argv[1:] if fn.endswith('.root')]

fns = ROOT.vector('string')()
for fname in fnames:
  fns.push_back(fname)


mm = re.search("\d{4,8}", fnames[0])
run = int(mm.group(0))


mpro = 0.938
E0 = 10.6
if run>=5674 and run<=5870:
  E0 = 7.54626
elif run>5870 and run<=6000:
  E0 = 6.53536


phistring = "double E0=" + str(E0) + ", Mpro=0.938;" + """

TLorentzVector beam(0,0,E0,E0);
TLorentzVector targ(0,0,0,Mpro);
TLorentzVector ele,pro;
ele.SetXYZM(ex,ey,ez,0);
pro.SetXYZM(px,py,pz,Mpro);

double phi = TMath::RadToDeg()*fabs(ele.Phi()-pro.Phi());

TLorentzVector wvec = beam+targ-ele;
double ww = wvec.M();

TVector3 zaxis(0,0,1);
TLorentzVector epx = wvec-pro;
double thgz = TMath::RadToDeg()*epx.Angle(zaxis);
double thge = TMath::RadToDeg()*epx.Angle(ele.Vect());

double e0 = 0.938*(-1+1/(tan(ele.Theta()/2)*tan(pro.Theta())));

double pe0 = (E0*Mpro)/(E0 + Mpro - E0*cos(ele.Theta()));
double dpe0 = pe0 - ele.P();

double pp0 = sqrt(E0 - pe0)*sqrt(E0 + 2*Mpro - pe0);
double dpp0 = pp0 - pro.P();

double pe1 = (e0*Mpro)/(e0 + Mpro - e0*cos(ele.Theta()));
double dpe1 = pe1 - ele.P();

double pp1 = sqrt(e0 - pe1)*sqrt(e0 + 2*Mpro - pe1);
double dpp1 = pp1 - pro.P();

double the = ele.Theta()*TMath::RadToDeg();
double thp = pro.Theta()*TMath::RadToDeg();
double thp0 = atan(Mpro/((E0 + Mpro)*tan(ele.Theta()/2.)))*TMath::RadToDeg();
double dthp0 = thp0-thp;

//double dpt = ele.Pt() - pro.Pt();
std::array<double,17> v{phi,ww,e0,ele.P(),dpe0,dpe1,pro.P(),dpp0,dpp1,thge,thgz,pe0,pe1,the,thp,thp0,dthp0};
return v;
"""

rdf = ROOT.RDataFrame("h22", fns)\
        .Filter("idet=="+str(idet))\
        .Define("vals",phistring)\
        .Define('fi','vals[0]')\
        .Define('ww','vals[1]')\
        .Define('e0','vals[2]')\
        .Define('pe','vals[3]')\
        .Define('dpe0','vals[4]')\
        .Define('dpe1','vals[5]')\
        .Define('pp','vals[6]')\
        .Define('dpp0','vals[7]')\
        .Define('dpp1','vals[8]')\
        .Define('thge','vals[9]')\
        .Define('thgz','vals[10]')\
        .Define('pe0','vals[11]')\
        .Define('pe1','vals[12]')\
        .Define('the','vals[13]')\
        .Define('thp','vals[14]')\
        .Define('thp0','vals[15]')\
        .Define('dthp0','vals[16]')

hfiw = [rdf.Filter('sec=='+str(i)).Histo2D(('hfiw'+str(i),'W vs #phi_{ep} (sec '+str(i)+');#phi_{ep} [#circ];W [GeV]',200,160,200,100,0,5),'fi','ww') for i in range(1,7)]
hfis = []

mus=[]
#c1 = ROOT.TCanvas("c1",'c1',1100,800)
#c1.Print('test.pdf[')
for isec in range(6):
  iw = hfiw[isec].GetYaxis().FindBin(1.5)
  h1 = hfiw[isec].ProjectionX("hfi"+str(isec+1), 0, iw)
  h1.SetTitle("#Delta#phi_{ep} (sec "+str(isec+1)+");#Delta#phi_{ep} [#circ]")
  hfis.append(h1)
  h1.Draw()
  f1 = ROOT.TF1("f1","gaus(0)+pol2(3)",177,183)
  f1.SetLineWidth(1)
  pars = []

  for isig in range(1,15):
    f1.SetParameters(h1.GetMaximum(),180,float(isig)/10,1,1,1)
    #print([f1.GetParameter(ipar) for ipar in range(6)])
    f1.SetParLimits(2,0,2)
    h1.Fit(f1,"RQ")
    if len(pars)==0 or pars[0]>f1.GetChisquare():
      pars = [f1.GetChisquare()]
      for ipar in range(6):
        pars.append(f1.GetParameter(ipar))
  for ipar in range(6):
    f1.SetParameter(ipar, pars[ipar+1])
  h1.Fit(f1,"RQ")
  #c1.Print('test.pdf')

  mu,sig = [f1.GetParameter(i) for i in range(1,3)]
  mus.extend([mu,sig])
#c1.Print('test.pdf]')

filterstr = "std::array<double, 12> phim{"+",".join(str(m) for m in mus)+"};return fabs(fi-phim[sec*2-2])<2*phim[sec*2-1];"

elasdf = rdf.Filter(filterstr)

cols = elasdf.GetColumnNames()
ivals = list(cols).index('vals')
cols.erase(cols.begin()+ivals)

fname = 'snap.{}.{:d}GeV.{:03d}.root'.format(detname, int(E0), len(fnames))
elasdf.Snapshot('h22',fname,cols)

ff = ROOT.TFile(fname,'update')
for i in range(6):
  hfiw[i].Write()
  hfis[i].Write()
ff.Close()
