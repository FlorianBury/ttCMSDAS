import os,sys
basepath = os.path.abspath(__file__).rsplit('/ttCMSDAS/',1)[0]+'/ttCMSDAS/'
sys.path.append(basepath)

from framework.analysis import analysis
from framework.functions import DeltaPhi, DiPt, InvMass, lepton, jet
from ROOT.TMath import Sqrt as sqrt
from ROOT import *

################ Analysis
class ttLeptonJet(analysis):
  def init(self):
    # Load SF files
    if not self.isData:
      self.LoadHisto('MuonIsoSF', basepath+'./inputs/MuonISO.root', 'NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta') # pt, abseta
      self.LoadHisto('MuonIdSF',  basepath+'./inputs/MuonID.root',  'NUM_TightID_DEN_genTracks_pt_abseta') # pt, abseta
      self.LoadHisto('ElecSF',    basepath+'./inputs/ElecTightCBid94X.root',  'EGamma_SF2D') # eta, pt

    # Objects for the analysis
    self.selLeptons = []
    self.selJets = []
    self.selBJets = []
    self.pmet = TLorentzVector()

    # Create output histograms
    self.categories = ["signal", "qcd", "signal_nomet", "qcd_nomet"]
    for selN in self.categories:
      self.myCreateTH1F(selN,"LepPt",   "", 24, 0, 120)
      self.myCreateTH1F(selN,"LepEta",  "", 50, -2.5, 2.5)
      self.myCreateTH1F(selN,"Bjet0_pt",  "", 30, 0, 150)
      self.myCreateTH1F(selN,"Bjet1_pt",  "", 30, 0, 150)
      self.myCreateTH1F(selN,"jet0_pt",  "", 30, 0, 150)
      self.myCreateTH1F(selN,"jet1_pt",  "", 30, 0, 150)
      self.myCreateTH1F(selN,"Bjet0_eta",  "", 50, -2.5, 2.5)
      self.myCreateTH1F(selN,"Bjet1_eta",  "", 50, -2.5, 2.5)
      self.myCreateTH1F(selN,"jet0_eta",  "", 50, -2.5, 2.5)
      self.myCreateTH1F(selN,"jet1_eta",  "", 50, -2.5, 2.5)
      self.myCreateTH1F(selN,"Bjets_InvMass",  "", 100, 0,500)
      self.myCreateTH1F(selN,"jets_InvMass",  "", 100, 0,500)
      self.myCreateTH1F(selN,"Bjets_DeltaPhi",  "", 100, -3.14/2,3.14/2)
      self.myCreateTH1F(selN,"Bjets_DeltaPt",  "", 50,0,200)
      self.myCreateTH1F(selN,"met_pt",  "", 50,0,200)

  def myCreateTH1F(self, sel, var, title, nBin, xMin, xMax):
    self.CreateTH1F("{0}_{1}".format(sel, var), title, nBin, xMin, xMax)
 
  def resetObjects(self):
    ''' Reset the list where the objects are stored '''
    self.selLeptons = []
    self.selJets = []
    self.selBJets = []
    self.pmet = TLorentzVector()

  def FillHistograms(self, sel, lepton, bjets, jets, pmet):
    ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
    if not len(lepton) >= 1: return # Just in case
    if not len(bjets) >= 2: return # Just in case
    if not len(jets) >= 2: return # Just in case
    self.weight = self.EventWeight * self.SFmuon * self.SFelec * self.PUSF

    # Re-calculate the observables
    lep_pt  = lepton.Pt()
    lep_eta = lepton.Eta()
    bjet0 = bjets[0] 
    bjet1 = bjets[1] 
    jet0 = jets[0] 
    jet1 = jets[1] 
    bjet0_pt = bjet0.Pt()
    bjet1_pt = bjet1.Pt()
    bjet0_eta = bjet0.Eta()
    bjet1_eta = bjet1.Eta()
    jet0_pt = jet0.Pt()
    jet1_pt = jet1.Pt()
    jet0_eta = jet0.Eta()
    jet1_eta = jet1.Eta()
    bjet_dphi  = DeltaPhi(bjet0, bjet1)
    mbb   = InvMass(bjet0, bjet1)
    mjj   = InvMass(jet0, jet1)
    bjet_dipt = DiPt(bjet0, bjet1)
    met_pt = pmet.Pt()
    
    ### Fill the histograms
    self.obj[sel+'_Lep_pt'].Fill(lep_pt, self.weight)
    self.obj[sel+'_Lep_eta'].Fill(lep_eta, self.weight)
    self.obj[sel+'_Bjet0_pt'].Fill(bjet0_pt, self.weight)
    self.obj[sel+'_Bjet0_eta'].Fill(bjet0_eta, self.weight)
    self.obj[sel+'_Bjet1_pt'].Fill(bjet1_pt, self.weight)
    self.obj[sel+'_Bjet1_eta'].Fill(bjet1_eta, self.weight)
    self.obj[sel+'_jet0_pt'].Fill(jet0_pt, self.weight)
    self.obj[sel+'_jet0_eta'].Fill(jet0_eta, self.weight)
    self.obj[sel+'_jet1_pt'].Fill(jet1_pt, self.weight)
    self.obj[sel+'_jet1_eta'].Fill(jet1_eta, self.weight)
    self.obj[sel+"_Bjets_InvMass"].Fill(mbb, self.weight)
    self.obj[sel+"_jets_InvMass"].Fill(mjj, self.weight)
    self.obj[sel+"_Bjets_DeltaPhi"].Fill(bjet_dphi, self.weight)
    self.obj[sel+"_Bjets_DeltaPt"].Fill(bjet_dipt, self.weight)
    self.obj[sel+"_met_pt"].Fill(met_pt, self.weight)

  def insideLoop(self, t):
    self.resetObjects()

    ### Lepton selection
    ###########################################
    if not self.isData: nGenLep = t.nGenDressedLepton 
    
    ##### Jets
    for i in range (t.nJet):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Jet_pt[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass[i])
      
      if p.Pt() < 30 or abs(p.Eta()) > 2.4:continue

      ## medium working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
      if t.Jet_btagDeepC > 0.4941 : 
        self.selBJets.append(p)
      else:
        self.selJets.append(p)

    ##### Muons
    QCDPass = False
    NoMetPass = False   
    for i in range(t.nMuon):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Muon_pt[i], t.Muon_eta[i], t.Muon_phi[i], t.Muon_mass[i])
      charge = t.Muon_charge[i]

      # Tight ID, tight ISO, RelIso04 < 0.15, tight IP
      if not t.Muon_tightId[i]: continue # Tight ID
      dxy = abs(t.Muon_dxy[i])
      dz  = abs(t.Muon_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 30 GeV, |eta| < 2.1
      if p.Pt() < 30 or abs(p.Eta()) > 2.1: continue

      if t.Muon_pfRelIso04_all[i] > 0.15 and t.MET_pt < 20 :
        self.selLeptons_QCD_noMet.append(lepton(p, charge, 13))
      elif t.Muon_pfRelIso04_all[i] < 0.15 and t.MET_pt < 20 :
        self.selLeptons__noMet.append(lepton(p, charge, 13))
      elif t.Muon_pfRelIso04_all[i] < 0.15 and t.MET_pt > 20 :
        self.selLeptons.append(lepton(p, charge, 13))
      else:
        self.selLeptons_QCD.append(lepton(p,charge,13))


    ##### Electrons
    for i in range(t.nElectron):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Electron_pt[i], t.Electron_eta[i], t.Electron_phi[i], t.Electron_mass[i])
      charge = t.Electron_charge[i]
      etaSC    = abs(p.Eta());
      convVeto = t.Electron_convVeto[i]

      # Tight cut-based Id, convVeto, RelIso03 tight, tight IP
      if not t.Electron_cutBased[i] >= 4: continue
      if not convVeto: continue
      relIso03 = t.Electron_pfRelIso03_all[i]
       
      dxy = abs(t.Electron_dxy[i])
      dz  = abs(t.Electron_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 30 GeV, |eta| < 2.1
      if p.Pt() < 30 or abs(p.Eta()) > 2.1: continue
      if   etaSC <= 1.479 and relIso03 > 0.0361: signalPass1 = True
      elif etaSC >  1.479 and relIso03 > 0.094:  signalPass2 = True
      if signalPass1 and signalPass2 and t.MET_pt < 20 :
        self.selLeptons_noMet.append(lepton(p, charge, 11))
      elif signalPass1 is False and signalPass2 is False and t.MET_pt < 20 :
        self.selLeptons_QCD_noMet.append(lepton(p, charge, 11))
      elif signalPass1 is False and signalPass2 is False and t.MET_pt > 20 :
        self.selLeptons_QCD.append(lepton(p, charge, 11))
      else:self.selLeptons.append(lepton(p, charge, 11))

    leps = self.selLeptons
    pts  = [lep.Pt() for lep in leps]
    self.selLeptons = [lep for _,lep in sorted(zip(pts,leps))]


    ##### MET 
    self.pmet = TLorentzVector()
    self.pmet.SetPtEtaPhiM(t.MET_pt, 0 , t.MET_phi, 0)

    ### Calculate the weights
    self.SFelec = 1; self.SFmuon = 1; self.SFelecErr = 0; self. SFmuonErr = 0
    if not self.isData:
      for lep in self.selLeptons:
        if lep.IsMuon():
          sf, err = self.GetSFandErr('MuonIsoSF, MuonIdSF', lep.Pt(), TMath.Abs(lep.Eta()))
          self.SFmuon*=sf
          self.SFmuonErr+=err*err
        else:
          sf, err = self.GetSFandErr('ElecSF', lep.Eta(), lep.Pt())
          self.SFelec*=sf
          self.SFelecErr+=err*err
      self.SFelecErr = sqrt(self.SFelecErr)
      self.SFmuonErr = sqrt(self.SFmuonErr)

    # PU SF --> PLEASE CHECK THAT THE WEIGHTS ARE IN THE TREES THAT YOU'RE USING!
    if not self.isData:
      self.PUSF   = t.puWeight
      self.PUUpSF = t.puWeightUp
      self.PUDoSF = t.puWeightDown
    else:
      self.PUSF   = 1; self.PUUpSF = 1; self.PUDoSF = 1

    ### Event selection
    ###########################################
    ### We need one lepton, two bjets and two light jets 
    ### Each of them must be at least 10 GeV
    if not len(leps) >= 1:      return 
    ### Fill the histograms
    print 'leptons ',self.selLeptons[:1]
    print 'bjets ',self.selBJets[:2]
    print 'jets ',self.selJets[:2]
    print 'met ',self.pmet
    

### SIGNAL
 
    ### One electron or muon and at least 4 jets
    ### veto events with more than one lepton.
    ### veto events containing leptons passing relaxed ID/ISO cuts. 
    if not len(leps) >= 1;      return
    if not len(selJets)>=2;     return
    if not len(selBJets)>=2;     return
     
    ### Fill the histograms
    self.FillHistograms("signal", self.selLeptons[:1], self.selBJets[:2], self.selJets[:2], self.pmet)


    ### QCD selection : non-prompt and less isolated leptons
    ### the control region : isolation cut is reversed 
    
    self.FillHistograms('qcd',self.selLeptons_QCD[:1], self.selBJets[:2], self.selJets[:2], self.pmet)
    self.FillHistograms('qcd_nomet',self.selLeptons_QCD_noMet[:1], self.selBJets[:2], self.selJets[:2], self.pmet)
    self.FillHistograms('signal_nomet',self.selLeptons_noMet[:1],self.selBJets[:2], self.selJets[:2], self.pmet)
