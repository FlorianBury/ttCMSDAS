from plotter.TopHistoReader import TopHistoReader, StackPlot, Process, WeightReader
from ROOT.TMath import Sqrt as sqrt
from ROOT import *

### Input and output
path = './temp/'
outpath = './outputs/'

### Definition of the processes
process = {
'VV'  : 'WZTo3LNU,WWTo2L2Nu',
'fake': 'WJetsToLNu,TTsemilep',
'tW'  : 'tW_noFullHad,  tbarW_noFullHad',
'DY'  : 'DYJetsToLL_M_10to50,DYJetsToLL_MLL50',
'tt'  : 'TT',
'data': 'HighEGJet,SingleMuon, DoubleMuon'}
prk = ['VV', 'fake', 'DY', 'tt', 'tW']

### Definition of colors for the processes
colors ={
'VV'  : kTeal+5,
'fake': kGray+2,
'tW'  : kOrange+1,
'DY'  : kAzure+2,
'tt'  : kRed+1,
'data': 1}

######################################################################################
### Plots
def DrawStack(out, *listOfPlots):
  ''' Draw some stack plots 
      Example:
      DrawStack(['NJets_ElMu_dilepton', 'Jet multiplicity'], ['DYMass_MuMu_dilepton', 'M_{#mu#mu}'])
  '''
  s = StackPlot(path)
  s.SetVerbose(1)
  s.SetOutPath(out)
  s.SetLumi(296.08)
  for pr in prk: s.AddProcess(pr, process[pr], colors[pr])
  s.AddData(process['data'])
  s.AddToSyst('ElecEffUp, ElecEffDown, MuonEffUp,MuonEffDown,PUUp, PUDown,JESUp,JESDown')
  s.SetRatioMin(0.5); s.SetRatioMax(1.5)
  for p in listOfPlots:
    if not isinstance(p, list): p = [p]
    s.DrawStack(p[0], p[1] if len(p) >= 2 else '', p[2] if len(p) >= 3 else '', p[3] if len(p) >= 4 else 1)

levels   = ['dilepton', 'MET','2jets']
channels = ['ElMu','ElEl', 'MuMu']
# Plots

for ch in channels:
  if   ch == 'MuMu': lepstr = '#mu#mu'
  elif ch == 'ElMu': lepstr = 'e#mu'
  elif ch == 'ElEl': lepstr = 'ee'
  for lev in levels:
    if lev == '1btag' and ch != 'ElMu': continue
    if ch != 'ElMu' and lev == 'dilepton': DrawStack(outpath+'leptons',[GetName('DYMass',ch,lev),  'M_{'+lepstr+'} (GeV)'])
    if lev == 'MET' and ch == 'ElMu': continue
    if lev == 'dilepton':
      DrawStack(outpath + 'global', 
      [GetName('NJets',ch,lev), 'Jet multiplicity'],
      [GetName('Vtx',ch,lev),  'Reconstructed primary verteces'])
    DrawStack(outpath + 'global',
      [GetName('HT',ch,lev),  'H_{T} (GeV)','Events / 40 GeV',8],
      [GetName('MET',ch,lev), 'Missing p_{T} (GeV)','Events / 15 GeV',3])
    DrawStack(outpath + 'leptons',
      [GetName('Lep0Pt',ch,lev),  'Leading lepton p_{T} (GeV)','Events / 10 GeV',2],
      [GetName('Lep1Pt',ch,lev),  'Subleading lepton p_{T} (GeV)','Events / 10 GeV',2],
      [GetName('Lep0Eta',ch,lev),  'Leading lepton #eta','',5],
      [GetName('Lep1Eta',ch,lev),  'Subleading lepton #eta','',5],
      [GetName('DilepPt',ch,lev),  'p_{T}^{'+lepstr+'} (GeV)','Events / 20 GeV',4],
      [GetName('InvMass',ch,lev),  'M_{'+lepstr+'} (GeV)','Events / 30 GeV',6],
      [GetName('DeltaPhi',ch,lev),  '#Delta#phi('+lepstr+')','',2])
    DrawStack(outpath + 'jets',
      [GetName('Jet0Pt',ch,lev),  'Leading Jet p_{T} (GeV)','Events / 30 GeV',6],
      [GetName('Jet1Pt',ch,lev),  'Subleading Jet p_{T} (GeV)','Evetns / 25 GeV',5],
      [GetName('Jet0Eta',ch,lev),  'Leading Jet #eta','',5],
      [GetName('Jet1Eta',ch,lev),  'Subleading Jet #eta','',5])
