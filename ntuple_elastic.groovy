#!/home/kenjo/.groovy/coatjava/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import org.jlab.jroot.ROOTFile
import my.Sugar


Sugar.enable()
/////////////////////////////////////////////////////////////////////

def hists = new ConcurrentHashMap()

def target = LorentzVector.withPID(2212,0,0,0)

def hdfi = {new H1F("$it","$it",400,150,210)}
def hw = {new H1F("$it","$it",200,0.5,4.5)}
def hmm2 = {new H1F("$it","$it",200,-1.5,1.5)}

def banknames = ['REC::Particle','REC::Calorimeter']

def ff = new ROOTFile("elastic_"+args[0].split('/').last().replace(".hipo", ".root"))
def tt = ff.makeNtuple('h22','elastic','ex:ey:ez:px:py:pz:idet:sec')

args.each{fname->
  def mm = fname.split('/')[-1] =~ /\d{4,7}/
  def run = mm[0].toInteger()

  def E0 = 10.6
  if(run>=5674 && run<=5870) E0 = 7.54626
  else if(run>5870 && run<=6000) E0 = 6.53536
  def beam = LorentzVector.withPID(11,0,0,E0)


  def reader = new HipoDataSource()
  reader.open(fname)

  while(reader.hasEvent()) {
    def event = reader.getNextEvent()

    if(banknames.every{event.hasBank(it)}) {
      def (partb,ecb) = banknames.collect{event.getBank(it)}

      (0..<partb.rows()).find{partb.getInt('pid',it)==11 && partb.getShort('status',it)<0}
      ?.with{iele->
        def ipros = (0..<partb.rows()).findAll{partb.getInt('pid',it)==2212}
        return ipros.size()==1 ? [iele,ipros[0]] : null
      }?.with{iele,ipro->
        def iec = (0..<ecb.rows()).find{ecb.getShort("pindex",it)==iele && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}
        return iec==null ? null : [iele,ipro,ecb.getByte('sector',iec)]
      }?.with{iele,ipro,esec->
        def idet = (partb.getShort('status',ipro)/1000).toInteger()
        def pdet = idet==2 ? 'FD' : 'CD'


        def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
        def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})

        def wvec = beam+target-ele
        def epx = beam+target-ele-pro

        def dfi = Math.toDegrees((ele.phi()-pro.phi()).abs())
        def ww = wvec.mass()

        [["",true],
         ["dfi.lt.2",(dfi-180).abs()<2],
         ["dfi.lt.2/w.lt.1.3",(dfi-180).abs()<2 && ww<1.3]
        ].findAll{it[1]}.each{dir,cut->
          hists.computeIfAbsent("$dir/hdfi:$pdet:$esec", hdfi).fill(dfi)
          hists.computeIfAbsent("$dir/hw:$pdet:$esec", hw).fill(ww)
          hists.computeIfAbsent("$dir/hmm2:$pdet:$esec", hmm2).fill(epx.mass2())
        }

        tt.fill(ele.px(), ele.py(), ele.pz(), pro.px(), pro.py(), pro.pz(), idet, esec)
      }
    }
  }

  reader.close()
}

tt.write()
hists.each{ff.writeDataSet(it.value)}
ff.close()

