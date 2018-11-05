// Macro to check that calibration is properly done: it weights events by the efficiency shape used to make the calibration to check if the bin distribution is flat
// Author : Javier Martin Blanco 10/04/2018

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <map>
//tree, hist, vector ...
#include <TROOT.h>
#include "TMath.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH1D.h>
#include <TF1.h>

//private setup
#include "../Utilities/treeUtils.h"

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
  {"step",{25.0,0.8156}},//0.631701 //0.509356 //0.8156 // Here the second parameter is the average efficiency in the inneficient region (given in the output of the centrality calibration)
  {"err",{13.439,0.811}},
  {"errPW",{11.36,13.5298,2.9839,0.2626}}
};


void checkCalibration(const char* inputFile = "../InputFiles/inputFiles.txt",
                      const char* effName = "step", // effName can be "step", "err" and "errPW"
                      const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1",
                      const char* sampleLabel = "Data_XeXe"
                      )
{
  // Read files and chain TTrees
  TChain* tree = readFiles(inputFile);
  //
  
  unsigned int nEvents = tree->GetEntries();
  if (nEvents>0) cout << "[INFO] Number of events to be processed = " << nEvents << endl;
  else {cout << "[ERROR] No events to be processed" << endl; return;}
  
  // Load needed branches
  if (!loadBranches(triggerName,tree)) {cout << "[ERROR] Error loading branches" << endl; return;}
  //
  
  TF1* fEff = effFfunc.at(effName);
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
  
  
  TH1::SetDefaultSumw2(true);
  
  int nBin = 200;
  double xMin = 0.;
  double xMax = 200;
  
  TH1D* hMB = new TH1D("hCentBin","", nBin,xMin,xMax);
  TH1D* hMB_w = new TH1D("hCentBin_weighted","", nBin,xMin,xMax);

  // Event loop
  for(unsigned int iev = 0; iev < nEvents; iev++)
  {
    if(iev%5000 == 0) cout<<"Processing event: " << iev << " / " << nEvents << endl;
    tree->GetEntry(iev);
    
    if ((HLT_Trigger_part1 || HLT_Trigger_part2 || HLT_Trigger_part3 || HLT_Trigger_part4 || HLT_Trigger_part5 ||
         HLT_Trigger_part6 || HLT_Trigger_part7 || HLT_Trigger_part8 || HLT_Trigger_part9 || HLT_Trigger_part10 ||
         HLT_Trigger_part11 || HLT_Trigger_part12 || HLT_Trigger_part13 || HLT_Trigger_part14 || HLT_Trigger_part15 ||
         HLT_Trigger_part16 || HLT_Trigger_part17 || HLT_Trigger_part18 || HLT_Trigger_part19 || HLT_Trigger_part20) && pprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter3)
    {
      // Lumi section selection
      if (run != 304899 && run != 304906) continue;
      if (run == 304899 && (lumi < 33 || lumi > 200)) continue;
      if (run == 304906 && (lumi < 1 || lumi > 731)) continue;
      
      double eff = 1./fEff->Eval(hf);
      if (eff > 1 && eff < 0.) eff = 1.;
      
      hMB->Fill(hiBin);
      hMB_w->Fill(hiBin,eff);
    }
  }
  
  // save results
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/Results";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  TFile *f = new TFile(Form("%s/centralityCheck_%s.root",saveDIR.c_str(),sampleLabel),"RECREATE");
  hMB->Write("hbin", TObject::kOverwrite | TObject::kSingleKey);
  hMB_w->Write("hbin_weighted", TObject::kOverwrite | TObject::kSingleKey);
  f->Close();
}
