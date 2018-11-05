#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

using namespace std;

bool descend(float i,float j) { return (i>j); }

double dEffErrF(double* x, double* par) // Use "err"
{
  // Error function
  
  return 0.5*(1+TMath::Erf((x[0]-par[0])/(TMath::Sqrt(x[0])*par[1])));
}

double dEffErrFPW(double* x, double* par) // Use "errPW"
{
  // Piecewise function changing at par[0]
  // Note that in this function par[2] and par[3] must be correlated in this way: par[3] = par[2]/par[0], to ensure the continuity
  
  if (x[0] <= par[0]) return 0.5*(1+TMath::Erf((x[0]-par[1])/par[2]));
  else return 0.5*(1+TMath::Erf((x[0]-par[1])/(par[3]*x[0])));
}

double dEffStepF(double* x, double* par) // Use "step"
{
  // Step function changing to 1
  // par[0] is the efficiency for x <= par[0]
  
  if (x[0] <= par[0]) return par[1];
  else return 1.0;
}

TF1* fEffStepF = new TF1("hfEfficiency",dEffStepF,0.,5000.,2);
TF1* fEffErrF = new TF1("hfEfficiency",dEffErrF,0.,5000.,2);
TF1* fEffErrFPW = new TF1("hfEfficiency",dEffErrFPW,0.,5000.,3);

const std::map< std::string , TF1* > effFfunc = {
  {"step",fEffStepF},
  {"err",fEffErrF},
  {"errPW",fEffErrFPW}
};

const std::map< std::string , std::vector<double> > effFpars = {
  {"step",{25.0,0.9476}},
  {"err",{13.439,0.811}},
  {"errPW",{11.36,13.5298,2.9839,0.2626}}
};

const std::map< std::string , std::string > effFstring = {
  {"step","if (x <= par[0]) eff = par[1] ; else eff= 1.0"},
  {"err","eff = 0.5*(1+TMath::Erf((x-par[0])/(TMath::Sqrt(x)*par[1])))"},
  {"errPW","if (x <= par[0]) eff = 0.5*(1+TMath::Erf((x-par[1])/par[2])) ; else eff = 0.5*(1+TMath::Erf((x-par[1])/(par[3]*x)))"}
};

void makeDataCentralityTable(int nbins = 200, const string label = "HFtowers",
                             const char * tag = "CentralityTable_HFtowers200_DataXeXe_effstep947a25_run2v941x02_offline",
                             bool useEffFunc = true, const char* effName = "step",double intEff = -1., // effName can be "step", "err" and "errPW"
                             bool geomInfo = true){
  
  TH1D::SetDefaultSumw2();
  
  TString inFileName = "files_XeXeData.txt";
  char line[1024];
  ifstream in(inFileName);
  TChain * t = new TChain("hiEvtAnalyzer/HiTree","");
  TChain * tskimanalysis = new TChain("skimanalysis/HltTree","");
  TChain * thltanalysis = new TChain("hltanalysisReco/HltTree","");
  while (in.getline(line,1024,'\n'))
  {
    t->Add(line);
    tskimanalysis->Add(line);
    thltanalysis->Add(line);
  }
  
  t->AddFriend(tskimanalysis);
  t->AddFriend(thltanalysis);
  
  TFile *outFile = new TFile("CentralityTable_HFtowers200_DataXeXe_effstep947a25_d20180412_v1.root","recreate");
  
  ofstream txtfile("output_DataXeXe_effstep947a25_d20180412_v1.txt");
  txtfile << "Input tree: " << inFileName << endl;
  txtfile << "Tag name: " << tag << endl;
  
  TDirectory *dir = outFile->mkdir(tag);
  dir->cd();
  TNtuple * nt = new TNtuple("nt","","value");
  
  const int runNum = 1;
  CentralityBins * bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);
  
  // Determine the binning to be used from input
  bool binHF = label.compare("HFtowers") == 0;
  bool binHFplus = label.compare("HFtowersPlus") == 0;
  bool binHFminus = label.compare("HFtowersMinus") == 0;
  bool binHFplusTrunc = label.compare("HFtowersPlusTrunc") == 0;
  bool binHFminusTrunc = label.compare("HFtowersMinusTrunc") == 0;
  bool binNpix = label.compare("PixelHits") == 0;
  bool binNpixTrks = label.compare("PixelTracks") == 0;
  bool binNtrks = label.compare("Tracks") == 0;
  
  // Deffine efficiency to weight events
  TF1* fEff(0x0);
  TParameter<double>* gEff(0x0);
  double effPrime = -1;
  if (useEffFunc)
  {
    fEff = effFfunc.at(effName);
    if (!fEff)
    {
      cout << "[ERROR] No efficiency function could be defined " << endl;
      return;
    }
    else cout << "[INFO] Using efficiency function ..." << endl;

    int parS = effFpars.at(effName).size();
    for (int i = 0 ; i < parS ; i++)
    {
      fEff->FixParameter(i,effFpars.at(effName).at(i));
    }
    
    if (!strcmp(effName,"step"))
    {
      cout << "Computing average efficiency for step function" << endl;
      
      const char* selection = "(HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1) && pPAprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter3";
      
      double ntot = (1./effFpars.at("step").at(1))*t->GetEntries(selection);
      double nlow = t->GetEntries(Form("%s && hiHF<=%f",selection,effFpars.at("step").at(0)));
      double nhigh = t->GetEntries(Form("%s && hiHF>%f",selection,effFpars.at("step").at(0)));
      effPrime = nlow / (ntot - nhigh);
      
      cout << "Average efficiency for  x <= " << effFpars.at("step").at(0) << " = " << effPrime << endl;
      
      fEff->FixParameter(1,effPrime);
    }
  }
  else
  {
    if (intEff>0 && intEff<=1) gEff = new TParameter<double>("eff",intEff);
    else
    {
      cout << "[ERROR] Input integrated efficiency is not correct, must be (0,1]" << fEff << endl;
      return;
    }
    cout << "Using global efficiency of " << gEff->GetVal() << endl;
  }
  
  
  //Here we need the default Glauber for 2.76 or 5 TeV
  TFile * inputMCfile(0x0);
  CentralityBins* inputMCtable(0x0);
  if (geomInfo)
  {
    inputMCfile= TFile::Open("/afs/cern.ch/work/j/jmartinb/private/Forest/CMSSW_9_4_1/src/HeavyIonsAnalysis/CentralityAnalysis/tools/CentralityTable_HFtowers200_XeXe5p44TeVEPOSLHC_d20180129_v1.root");
    inputMCtable = (CentralityBins*)inputMCfile->Get("CentralityTable_HFtowers200_XeXe5p44TeVEPOSLHC_v941x01_mc/run1");
  }
  
  
  double binboundaries[nbins+1];
  vector<float> values;
  
  float hf, hfplus, hfpluseta4, hfminuseta4, hfminus, hfhit, ee, eb, zdc, zdcplus, zdcminus;
  int HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1,HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1, HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1;
  int run, lumi, npix, npixtrks, ntrks ,pPAprimaryVertexFilter, pBeamScrapingFilter, phfCoincFilter3;
  t->SetBranchAddress("run",    &run);
  t->SetBranchAddress("lumi",   &lumi);
  t->SetBranchAddress("hiHF",           &hf);
  t->SetBranchAddress("hiHFplus",       &hfplus);
  t->SetBranchAddress("hiHFplusEta4",   &hfpluseta4);
  t->SetBranchAddress("hiHFminus",      &hfminus);
  t->SetBranchAddress("hiHFminusEta4",  &hfminuseta4);
  t->SetBranchAddress("hiHFhit",        &hfhit);
  t->SetBranchAddress("hiZDC",          &zdc);
  t->SetBranchAddress("hiZDCplus",      &zdcplus);
  t->SetBranchAddress("hiZDCminus",     &zdcminus);
  t->SetBranchAddress("hiEE",           &ee);
  t->SetBranchAddress("hiEB",           &eb);
  t->SetBranchAddress("hiNpix",         &npix);
  t->SetBranchAddress("hiNpixelTracks", &npixtrks);
  t->SetBranchAddress("hiNtracks",      &ntrks);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1);
  t->SetBranchAddress("HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1", &HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1);
  t->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter);
  t->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);
  t->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3);
  
  unsigned int Nevents = t->GetEntries();
  
  double totalXsec(0);
  int passed(0);
  for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if(iev%5000 == 0) cout<<"Processing event: " << iev << " / " << Nevents << endl;
    t->GetEntry(iev);
    
    if ((HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1 || HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1) && pPAprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter3 )
    {
      passed++;
      
      float parameter = -1;
      if(binHF) parameter = hf;
      if(binHFplus) parameter = hfplus;
      if(binHFminus) parameter = hfminus;
      if(binHFplusTrunc) parameter = hfpluseta4;
      if(binHFminusTrunc) parameter = hfminuseta4;
      if(binNpix) parameter = npix;
      if(binNpixTrks) parameter = npixtrks;
      if(binNtrks) parameter = ntrks;
      
      values.push_back(parameter);
      nt->Fill(parameter);
      
      if (useEffFunc)
      {
        double eff = fEff->Eval(parameter);
//        cout << "eff(" << parameter << ") = " << eff << endl;
        if (eff <= 1 && eff > 0.) totalXsec += 1./eff;
        else totalXsec += 1.;
      }
    }
  }
  
  cout << endl;
  cout << "Selected events = " << passed << endl;
  cout << "Selected weighed events = " << totalXsec << endl;
  cout << "Total efficiency = " << (double)(passed/totalXsec) << endl;
  
  sort(values.begin(),values.end());
  
  txtfile << "Number of events = " << passed << endl;
  txtfile << "Number of weighted events = " << totalXsec << endl;
  txtfile << "Total efficiency = " << (double)(passed/totalXsec) << endl;
  txtfile << endl;
  txtfile << "-------------------------------------" << endl;
  if (useEffFunc)
  {
    txtfile << "EFficiency shape is " << effName << " :" << endl;
    txtfile << effFstring.at(effName) << endl;
    txtfile << "with pars : " << endl;
    for (int i = 0 ; i < effFpars.at(effName).size() ; i++)
    {
      txtfile << Form("par[%d] = %f",i,fEff->GetParameter(i)) << endl;
    }
    
    if (effPrime>0) txtfile << "Average efficiency for  x <= " << effFpars.at("step").at(0) << " = " << effPrime << endl;
  }
  txtfile << endl;
  txtfile << "-------------------------------------" << endl;
  txtfile << label.data() << " based cuts are: " << endl;
  txtfile << "(";
  
  unsigned int size = values.size();
  binboundaries[nbins] = values[size-1];
  cout << "Events per bin = " << totalXsec/(double)(nbins) << endl;
  if (useEffFunc)
  {
    binboundaries[0] = 0.;
    
    double integral = 0;
    int currentbin = 1;
    for(unsigned int iev = 0; iev < size; iev++)
    {
      double val = values[iev];
      double eff = fEff->Eval(val);
      if (eff <= 1 && eff > 0.) integral += 1. / eff;
      else integral += 1.;
      
      if(integral > ((double)(currentbin))*(totalXsec/(double)(nbins)))
      {
        cout << "current bin = " << currentbin << " ; integral = " << integral << " ; sum = " << ((double)(currentbin))*(totalXsec/(double)(nbins)) << endl;
        binboundaries[currentbin] = val;
        currentbin++;
      }
    }
    cout << currentbin << endl;
  }
  else // This way assumes all inefficiency is in the most peripheral bins
  {
    double EFF = gEff->GetVal();
    for(int i = 0; i < nbins; i++) {
      int entry = (int)( (float)(i)*((float)(size)/EFF/(float)(nbins)) - (float)(size)*(1. - EFF)/EFF );
      if(entry < 0 || i == 0) binboundaries[i] = 0;
      else binboundaries[i] = values[entry];
      if(binboundaries[i] < 0) { binboundaries[i] = 0; cout << "*"; }
    }
  }
  
  for(int i = 0; i < nbins; i++) {
    if(binboundaries[i] < 0) binboundaries[i] = 0;
    txtfile << binboundaries[i] << ", ";
  }
  txtfile << binboundaries[nbins] << ")" << endl;
  txtfile << endl;
  
  txtfile<<"-------------------------------------"<<endl;
  txtfile<<"# Bin NpartMean NpartSigma NcollMean NcollSigma bMean bSigma BinEdge"<<endl;
  for(int i = 0; i < nbins; i++){
    int ii = nbins-i;
    
    if (inputMCtable)
    {
      bins->table_[i].n_part_mean = inputMCtable->NpartMeanOfBin(i);
      bins->table_[i].n_part_var = inputMCtable->NpartSigmaOfBin(i);
      bins->table_[i].n_coll_mean = inputMCtable->NcollMeanOfBin(i);
      bins->table_[i].n_coll_var = inputMCtable->NcollSigmaOfBin(i);
      bins->table_[i].b_mean = inputMCtable->bMeanOfBin(i);
      bins->table_[i].b_var = inputMCtable->bSigmaOfBin(i);
      bins->table_[i].n_hard_mean = inputMCtable->NhardMeanOfBin(i);
      bins->table_[i].n_hard_var = inputMCtable->NhardSigmaOfBin(i);
      bins->table_[i].ecc2_mean  = inputMCtable->eccentricityMeanOfBin(i);
      bins->table_[i].ecc2_var = inputMCtable->eccentricitySigmaOfBin(i);
    }
    else
    {
      bins->table_[i].n_part_mean = -99;
      bins->table_[i].n_part_var = -99;
      bins->table_[i].n_coll_mean = -99;
      bins->table_[i].n_coll_var = -99;
      bins->table_[i].b_mean = -99;
      bins->table_[i].b_var = -99;
      bins->table_[i].n_hard_mean = -99;
      bins->table_[i].n_hard_var = -99;
      bins->table_[i].ecc2_mean  = -99;
      bins->table_[i].ecc2_var = -99;
    }
    
    bins->table_[i].bin_edge = binboundaries[ii-1];
    
    txtfile << i << " " << bins->table_[i].n_part_mean << " " << bins->table_[i].n_part_var << " " << bins->table_[i].n_coll_mean << " " << bins->table_[i].n_coll_var << " " <<bins->table_[i].b_mean << " " << bins->table_[i].b_var << " " << bins->table_[i].n_hard_mean << " " << bins->table_[i].n_hard_var << " " << bins->table_[i].bin_edge << " " << endl;
  }
  txtfile << endl;
  txtfile<<"-------------------------------------"<<endl;
  
  outFile->cd();
  dir->cd();
  bins->Write();
  nt->Write();
//  bins->Delete();
  outFile->Write();
  outFile->Close();
  txtfile.close();
  
}
