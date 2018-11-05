// comparison centrality objects between DATA 2016 and DATA 2015
// Author : Javier Martin Blanco 10/11/2016

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <cstdlib>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>

//private setup
#include "../Utilities/treeUtils.h"

void mbDistributions(const char* inputFile = "../InputFiles/inputFiles.txt",
                    const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1",
                    const char* sampleLabel = "STARlight_single"
                    )
{
  // Read files and chain TTrees
  TChain* tree = readFiles(inputFile);
  //
  
  unsigned int nEvents = tree->GetEntries();
  if (nEvents>0) cout << "[INFO] Number of events to be processed = " << nEvents << endl;
  else {cout << "[ERROR] No events to be processed" << endl; return;}
  
  if (nEvents > 200000) nEvents = 200000;
    
  // Load needed branches
  if (!loadBranches(triggerName,tree)) {cout << "[ERROR] Error loading branches" << endl; return;}
  //
  
  
  // Definition of variable ranges and number of bins for plots
  double hiHFMax = 4000;
  double hiEMax = 3000;
  double hiEtMax = 500;
  double hiHFSSMax = 2000;
  double hiHFSSTruncMax = 80;
  double hiBinMax = 100;
  double hiHFhitMax = 10000;
  double hiETMax = 50;
  double hiEEMax = 50;
  double hiMBMax = 100;
  double hiNpixMax = 40000;
  double hiNpixelTracksMax = 500;
  double hiNtracksMax = 4000;
  double hiNtracksCutEtaMax = 350;
  double hiNtracksCutMax = 150;
  double hiZDCMax = 40000;
  
  int nBin = 100;
  
  const int nvar = 15;//26;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFminus", "hiHFasym", "leadingTowerE", "leadingTowerEp","leadingTowerEm", "leadingTowerET", "leadingTowerETp", "leadingTowerETm","leadingTowerEAsym","leadingTowerETAsym","hiNtracks","trkEta","trkPt"};//,"nTrk", "hiHFplusEta4","hiHFminusEta4","hiHFhit","hiNpix","zVtx","nVtx","trkDxy1","trkDz1","nTrkVtx","highPurity","normChi2Vtx","vtxDist2D"};
  
  TString tLabel(triggerName);
  if (tLabel.Contains("MinimumBias")) tLabel = "MB";
  else if (tLabel.Contains("ZeroBias")) tLabel = "ZB";
  else if (tLabel.Contains("UnpairedBunch")) tLabel = "UB";
  else if (tLabel.Contains("NotBptxOR")) tLabel = "NBOR";
  else {cout << "[ERROR] Unkown trigger label" << endl; return;}

    // Definition of histos to store distributions
  TObjArray* avarName = new TObjArray();
  for (int i = 0 ; i < nvar ; i++)
  {
    avarName->Add(new TObjString(varname[i]));
  }
  
  TObjArray* aSave = new TObjArray();
  
  TH1::SetDefaultSumw2(true);
  for (int i = 0 ; i < nvar ; i++)
  {
    TObjArray* arr = new TObjArray();
    arr->SetName(Form("histos_%s",varname[i]));
    
    double xMin = 0.;
    double xMax = 300;
    
    if( !strcmp(varname[i],"hiNtracks") ) {xMax = hiNtracksMax; nBin = (int)hiNtracksMax;}
    if( !strcmp(varname[i],"hiHFasym") || !strcmp(varname[i],"leadingTowerEAsym") || !strcmp(varname[i],"leadingTowerETAsym")) {xMin = -1.2; xMax = 1.2; nBin = 240;}
    if( !strcmp(varname[i],"hiHF")) {xMax = hiHFMax; nBin = (int)(hiHFMax);}
    if( !strcmp(varname[i],"leadingTowerE") || !strcmp(varname[i],"leadingTowerEp") || !strcmp(varname[i],"leadingTowerEm")) {xMax = hiEMax; nBin = (int)(hiEMax);}
    if( !strcmp(varname[i],"leadingTowerET") || !strcmp(varname[i],"leadingTowerETm") || !strcmp(varname[i],"leadingTowerETp")) {xMax = hiEtMax; nBin = (int)(hiEtMax);}
    if( !strcmp(varname[i],"hiHFplus") || !strcmp(varname[i],"hiHFminus")) {xMax = hiHFSSMax; nBin = (int)(hiHFSSMax/2.);}
    if( !strcmp(varname[i],"hiHFplusEta4") || !strcmp(varname[i],"hiHFminusEta4")) {xMax = hiHFSSTruncMax; nBin = (int)hiHFSSTruncMax*8;}
    if( !strcmp(varname[i],"hiHFhit") ){xMax = hiHFhitMax; nBin = (int)(hiHFhitMax/10.);}
    if( !strcmp(varname[i],"hiNpix") ){xMax = hiNpixMax; nBin = (int)(hiNpixMax/10.);}
    if( !strcmp(varname[i],"zVtx") ){xMin = -30.; xMax = 30.; nBin = 1200;}
    if( !strcmp(varname[i],"trkDxy1") ){xMin = -150; xMax = 150.; nBin = 1200;}
    if( !strcmp(varname[i],"trkDz1") ){xMin = -600; xMax = 600.; nBin = 4800;}
    if( !strcmp(varname[i],"nVtx") ){xMax = 20; nBin = 20;}
    if( !strcmp(varname[i],"nTrkVtx") ){xMax = 500; nBin = 500;}
    if( !strcmp(varname[i],"highPurity") ){xMax = 2; nBin = 2;}
    if( !strcmp(varname[i],"nTrk") ){xMax = 600; nBin = 600;}
    if( !strcmp(varname[i],"trkEta") ){xMin = -4. ; xMax = 4.; nBin = 40;}
    if( !strcmp(varname[i],"trkPt") ){xMin = 0. ; xMax = 1000.; nBin = 500;}
    if( !strcmp(varname[i],"normChi2Vtx") ){xMax = 200; nBin = 2000;}
    if( !strcmp(varname[i],"vtxDist2D") ){xMin = -3;xMax = 3; nBin = 1920;}

    
    TH1D* h_MB(0x0);
    if( !strcmp(varname[i],"trkPt") )
    {
      const int nBinTrk = 1000;
      Double_t xBin[nBinTrk+1];
      for (int j = 0 ; j < nBinTrk ; j++)
      {
        if (j == 0) xBin[j] = 0;
        if (j <= 30) xBin[j] = 0.1*j; //bw = 0.1
        if (j > 30 && j <= 40) xBin[j] = 0.1*j + 0.1*(j-30);  //bw = 0.2
        if (j > 40 && j <= 60) xBin[j] = 0.1*j + 0.1*(j-30) + 0.3*(j-40);  //bw = 0.5
        if (j > 60 && j <= 80) xBin[j] = 0.1*j + 0.1*(j-30) + 0.3*(j-40) + 0.5*(j-60);  //bw = 1
        if (j > 80) xBin[j] = 0.1*j + 0.1*(j-30) + 0.3*(j-40) + 0.5*(j-60) + (j-80);  //bw = 2
      }
      xBin[nBinTrk] = 0.1*nBinTrk + 0.1*(nBinTrk-30) + 0.3*(nBinTrk-40) + 0.5*(nBinTrk-60) + (nBinTrk-80);
      h_MB = new TH1D(Form("h_%s_%s",tLabel.Data(),varname[i]), Form(";%s;",varname[i]), nBinTrk,xBin);
    }
    else h_MB = new TH1D(Form("h_%s_%s",tLabel.Data(),varname[i]), Form(";%s;",varname[i]), nBin,xMin,xMax);
    TH1D* h = (TH1D*)h_MB->Clone(Form("h_%s",varname[i]));
    TH1D* h_BS = (TH1D*)h_MB->Clone(Form("h_BS_%s",varname[i]));
    TH1D* h_PV = (TH1D*)h_MB->Clone(Form("h_PV_%s",varname[i]));
    TH1D* h_HF1 = (TH1D*)h_MB->Clone(Form("h_HF1_%s",varname[i]));
    TH1D* h_HF2 = (TH1D*)h_MB->Clone(Form("h_HF2_%s",varname[i]));
    TH1D* h_HF3 = (TH1D*)h_MB->Clone(Form("h_HF3_%s",varname[i]));
    TH1D* h_HF4 = (TH1D*)h_MB->Clone(Form("h_HF4_%s",varname[i]));
    TH1D* h_HF5 = (TH1D*)h_MB->Clone(Form("h_HF5_%s",varname[i]));
    TH1D* h_BS_PV = (TH1D*)h_MB->Clone(Form("h_BS_PV_%s",varname[i]));
    TH1D* h_BS_PV_HF1 = (TH1D*)h_MB->Clone(Form("h_BS_PV_HF1_%s",varname[i]));
    TH1D* h_BS_PV_HF2 = (TH1D*)h_MB->Clone(Form("h_BS_PV_HF2_%s",varname[i]));
    TH1D* h_BS_PV_HF3 = (TH1D*)h_MB->Clone(Form("h_BS_PV_HF3_%s",varname[i]));
    TH1D* h_BS_PV_HF4 = (TH1D*)h_MB->Clone(Form("h_BS_PV_HF4_%s",varname[i]));
    TH1D* h_BS_PV_HF5 = (TH1D*)h_MB->Clone(Form("h_BS_PV_HF5_%s",varname[i]));
    
    TH1D* h_MB_BS = (TH1D*)h_MB->Clone(Form("h_%s_BS_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_PV = (TH1D*)h_MB->Clone(Form("h_%s_PV_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_HF1 = (TH1D*)h_MB->Clone(Form("h_%s_HF1_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_HF2 = (TH1D*)h_MB->Clone(Form("h_%s_HF2_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_HF3 = (TH1D*)h_MB->Clone(Form("h_%s_HF3_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_HF4 = (TH1D*)h_MB->Clone(Form("h_%s_HF4_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_PV = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_PV_HF1 = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_HF1_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_PV_HF2 = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_HF2_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_PV_HF3 = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_HF3_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_PV_HF4 = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_HF4_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_PV_HF5 = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_HF5_%s",tLabel.Data(),varname[i]));
    
    TH1D* h_MB_BS_PV_NotHF3 = (TH1D*)h_MB->Clone(Form("h_%s_BS_PV_NotHF3_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_NotPV = (TH1D*)h_MB->Clone(Form("h_%s_BS_NotPV_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_NotPV_HF3 = (TH1D*)h_MB->Clone(Form("h_%s_BS_NotPV_HF3_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_BS_NotPV_NotHF3 = (TH1D*)h_MB->Clone(Form("h_%s_BS_NotPV_NotHF3_%s",tLabel.Data(),varname[i]));
    
    TH1D* h_MB_NotBS = (TH1D*)h_MB->Clone(Form("h_%s_NotBS_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_NotPV = (TH1D*)h_MB->Clone(Form("h_%s_NotPV_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_NotHF1 = (TH1D*)h_MB->Clone(Form("h_%s_NotHF1_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_NotHF2 = (TH1D*)h_MB->Clone(Form("h_%s_NotHF2_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_NotHF3 = (TH1D*)h_MB->Clone(Form("h_%s_NotHF3_%s",tLabel.Data(),varname[i]));
    TH1D* h_MB_NotHF4 = (TH1D*)h_MB->Clone(Form("h_%s_NotHF4_%s",tLabel.Data(),varname[i]));
    
    TH2D* h_asym = new TH2D(Form("h_hiHFasym_vs_%s",varname[i]), Form(";%s;hiHFasym",varname[i]), nBin,xMin,xMax, 240, -1.2, 1.2);
    TH2D* h_asym_MB = (TH2D*)h_asym->Clone(Form("h_%s_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asym_MB_BS = (TH2D*)h_asym->Clone(Form("h_%s_BS_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
    TH2D* h_asym_MB_BS_PV = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asym_MB_BS_PV_HF1 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF1_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asym_MB_BS_PV_HF2 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF2_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
    TH2D* h_asym_MB_BS_PV_HF3 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF3_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asym_MB_BS_PV_HF4 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF4_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asym_MB_BS_PV_HF5 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF5_hiHFasym_vs_%s",tLabel.Data(),varname[i]));
    
    TH2D* h_asymLE = (TH2D*)h_asym->Clone(Form("h_hiHFasymLE_vs_%s",varname[i]));
    TH2D* h_asymLE_MB = (TH2D*)h_asym->Clone(Form("h_%s_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLE_MB_BS = (TH2D*)h_asym->Clone(Form("h_%s_BS_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
    TH2D* h_asymLE_MB_BS_PV = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLE_MB_BS_PV_HF1 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF1_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLE_MB_BS_PV_HF2 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF2_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
    TH2D* h_asymLE_MB_BS_PV_HF3 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF3_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLE_MB_BS_PV_HF4 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF4_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLE_MB_BS_PV_HF5 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF5_hiHFasymLE_vs_%s",tLabel.Data(),varname[i]));
    
    TH2D* h_asymLET = (TH2D*)h_asym->Clone(Form("h_hiHFasymLET_vs_%s",varname[i]));
    TH2D* h_asymLET_MB = (TH2D*)h_asym->Clone(Form("h_%s_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLET_MB_BS = (TH2D*)h_asym->Clone(Form("h_%s_BS_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
    TH2D* h_asymLET_MB_BS_PV = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLET_MB_BS_PV_HF1 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF1_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLET_MB_BS_PV_HF2 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF2_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
    TH2D* h_asymLET_MB_BS_PV_HF3 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF3_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLET_MB_BS_PV_HF4 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF4_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
//    TH2D* h_asymLET_MB_BS_PV_HF5 = (TH2D*)h_asym->Clone(Form("h_%s_BS_PV_HF5_hiHFasymLET_vs_%s",tLabel.Data(),varname[i]));
    
    arr->Add(h);
    arr->Add(h_BS);
    arr->Add(h_PV);
    arr->Add(h_HF1);
    arr->Add(h_HF2);
    arr->Add(h_HF3);
    arr->Add(h_HF4);
    arr->Add(h_HF5);
    arr->Add(h_BS_PV);
    arr->Add(h_BS_PV_HF1);
    arr->Add(h_BS_PV_HF2);
    arr->Add(h_BS_PV_HF3);
    arr->Add(h_BS_PV_HF4);
    arr->Add(h_BS_PV_HF5);
    
    arr->Add(h_MB);
    arr->Add(h_MB_BS);
    arr->Add(h_MB_PV);
    arr->Add(h_MB_HF1);
    arr->Add(h_MB_HF2);
    arr->Add(h_MB_HF3);
    arr->Add(h_MB_HF4);
    arr->Add(h_MB_BS_PV_NotHF3);
    arr->Add(h_MB_BS_NotPV);
    arr->Add(h_MB_BS_NotPV_HF3);
    arr->Add(h_MB_BS_NotPV_NotHF3);
    arr->Add(h_MB_BS_PV);
    arr->Add(h_MB_BS_PV_HF1);
    arr->Add(h_MB_BS_PV_HF2);
    arr->Add(h_MB_BS_PV_HF3);
    arr->Add(h_MB_BS_PV_HF4);
    arr->Add(h_MB_BS_PV_HF5);
    
    arr->Add(h_MB_NotBS);
    arr->Add(h_MB_NotPV);
    arr->Add(h_MB_NotHF1);
    arr->Add(h_MB_NotHF2);
    arr->Add(h_MB_NotHF3);
    arr->Add(h_MB_NotHF4);
    
    arr->Add(h_asym);
    arr->Add(h_asym_MB);
//    arr->Add(h_asym_MB_BS);
    arr->Add(h_asym_MB_BS_PV);
//    arr->Add(h_asym_MB_BS_PV_HF1);
//    arr->Add(h_asym_MB_BS_PV_HF2);
    arr->Add(h_asym_MB_BS_PV_HF3);
//    arr->Add(h_asym_MB_BS_PV_HF4);
//    arr->Add(h_asym_MB_BS_PV_HF5);
    
    arr->Add(h_asymLE);
    arr->Add(h_asymLE_MB);
//    arr->Add(h_asymLE_MB_BS);
    arr->Add(h_asymLE_MB_BS_PV);
//    arr->Add(h_asymLE_MB_BS_PV_HF1);
//    arr->Add(h_asymLE_MB_BS_PV_HF2);
    arr->Add(h_asymLE_MB_BS_PV_HF3);
//    arr->Add(h_asymLE_MB_BS_PV_HF4);
//    arr->Add(h_asymLE_MB_BS_PV_HF5);
    
    arr->Add(h_asymLET);
    arr->Add(h_asymLET_MB);
//    arr->Add(h_asymLET_MB_BS);
    arr->Add(h_asymLET_MB_BS_PV);
//    arr->Add(h_asymLET_MB_BS_PV_HF1);
//    arr->Add(h_asymLET_MB_BS_PV_HF2);
    arr->Add(h_asymLET_MB_BS_PV_HF3);
//    arr->Add(h_asymLET_MB_BS_PV_HF4);
//    arr->Add(h_asymLET_MB_BS_PV_HF5);

  
    aSave->Add(arr);
  }
  
  TObjArray* adummy(0x0);
  // Event loop
  for(unsigned int iev = 0; iev < nEvents; iev++)
  {
    if(iev%5000 == 0) cout<<"Processing event: " << iev << " / " << nEvents << endl;
    tree->GetEntry(iev);
  
    if (ProcessID == 102) continue; // Reject elastic events in MC
    
    // Compute Energy of leading towers for this event
    double leadingEp(0.),leadingETp(0.),leadingEm(0.),leadingETm(0.),leadingEAsym(0.),leadingETAsym(0.);
    
    vector<float> Ep;
    vector<float> Em;
    vector<float> ETp;
    vector<float> ETm;
    for (int tower = 0 ; tower < nTowers ; tower++)
    {
      bool isHF = (abs(ieta[tower])>29);
      if (isHF && eta[tower]>0)
      {
        Ep.push_back(towerE[tower]);
        ETp.push_back(towerET[tower]);
      }
      if (isHF && eta[tower]<0)
      {
        Em.push_back(towerE[tower]);
        ETm.push_back(towerET[tower]);
      }
    }
    
    sort(Ep.begin(),Ep.end());
    sort(Em.begin(),Em.end());
    sort(ETp.begin(),ETp.end());
    sort(ETm.begin(),ETm.end());
    
    leadingEp = Ep.at(Ep.size() -1);
    leadingEm = Em.at(Em.size() -1);
    leadingETp = ETp.at(ETp.size() -1);
    leadingETm = ETm.at(ETm.size() -1);
    
    leadingEAsym = (leadingEp - leadingEm)/(leadingEp + leadingEm);
    leadingETAsym = (leadingETp - leadingETm)/(leadingETp + leadingETm);
    //=================================
    
    for (int i = 0 ; i < nvar ; i++)
    {
      TString varN = static_cast<TObjString*>(avarName->At(i))->GetString();
      adummy = static_cast<TObjArray*>(aSave->FindObject(Form("histos_%s",varN.Data())));
      
      float parameter = -1;
      int loopPar = 1;
      if(!strcmp(varN.Data(),"hiHF")) parameter = hf;
      if(!strcmp(varN.Data(),"leadingTowerE")) parameter = leadingEp+leadingEm;
      if(!strcmp(varN.Data(),"leadingTowerET")) parameter = leadingETp+leadingETm;
      if(!strcmp(varN.Data(),"leadingTowerEp")) parameter = leadingEp;
      if(!strcmp(varN.Data(),"leadingTowerEm")) parameter = leadingEm;
      if(!strcmp(varN.Data(),"leadingTowerETp")) parameter = leadingETp;
      if(!strcmp(varN.Data(),"leadingTowerETm")) parameter = leadingETm;
      if(!strcmp(varN.Data(),"leadingTowerEAsym")) parameter = leadingEAsym;
      if(!strcmp(varN.Data(),"leadingTowerETAsym")) parameter = leadingETAsym;
      if(!strcmp(varN.Data(),"hiHFasym")) parameter = (hfplus-hfminus)/hf;
      if(!strcmp(varN.Data(),"hiHFplus")) parameter = hfplus;
      if(!strcmp(varN.Data(),"hiHFminus")) parameter = hfminus;
      if(!strcmp(varN.Data(),"hiHFplusEta4")) parameter = hfpluseta4;
      if(!strcmp(varN.Data(),"hiHFminusEta4")) parameter = hfminuseta4;
      if(!strcmp(varN.Data(),"hiNtracks")) parameter = hiNtrks;
      if(!strcmp(varN.Data(),"hiHFhit")) parameter = hfhit;
      if(!strcmp(varN.Data(),"hiNpix")) parameter = npix;
      if(!strcmp(varN.Data(),"nVtx")) parameter = nvtx;
      if(!strcmp(varN.Data(),"nTrk")) parameter = ntrk;
      if(!strcmp(varN.Data(),"zVtx")) parameter = zvtx[maxptvtx];
      if(!strcmp(varN.Data(),"trkDxy1")) parameter = trkdxy1[maxptvtx];
      if(!strcmp(varN.Data(),"trkDz1")) parameter = trkdz1[maxptvtx];
      if(!strcmp(varN.Data(),"nTrkVtx")) parameter = ntrkvtx[maxptvtx];
      if(!strcmp(varN.Data(),"highPurity")) parameter = (int)highpurity[maxptvtx];
      if(!strcmp(varN.Data(),"normChi2Vtx") )parameter = normchi2vtx[maxptvtx];
      if(!strcmp(varN.Data(),"vtxDist2D") )parameter = vtxdist2d[maxptvtx];
      if(!strcmp(varN.Data(),"trkEta") ) {loopPar = ntrk;}
      if(!strcmp(varN.Data(),"trkPt") ) {loopPar = ntrk;}
  
      for (int n = 0 ; n < loopPar ; n++ )
      {
        if(!strcmp(varN.Data(),"trkEta") ) parameter = trkEta[n];
        if(!strcmp(varN.Data(),"trkPt") )  parameter = trkPt[n];
        
        static_cast<TH1*>(adummy->FindObject(Form("h_%s",varN.Data())))->Fill(parameter); // Histogram without any selection
        static_cast<TH2*>(adummy->FindObject(Form("h_hiHFasym_vs_%s",varN.Data())))->Fill(parameter,(hfplus-hfminus)/hf); // Histogram without any selection
        static_cast<TH2*>(adummy->FindObject(Form("h_hiHFasymLE_vs_%s",varN.Data())))->Fill(parameter,leadingEAsym); // Histogram without any selection
        static_cast<TH2*>(adummy->FindObject(Form("h_hiHFasymLET_vs_%s",varN.Data())))->Fill(parameter,leadingETAsym); // Histogram without any selection
        
        if (pBeamScrapingFilter)
        {
          static_cast<TH1*>(adummy->FindObject(Form("h_BS_%s",varname[i])))->Fill(parameter);
          
          if (pprimaryVertexFilter)
          {
           static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_%s",varname[i])))->Fill(parameter);
            
            if (phfCoincFilter1) static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF1_%s",varname[i])))->Fill(parameter);
            if (phfCoincFilter2) static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF2_%s",varname[i])))->Fill(parameter);
            if (phfCoincFilter3) static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF3_%s",varname[i])))->Fill(parameter);
            if (phfCoincFilter4) static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF4_%s",varname[i])))->Fill(parameter);
            if (phfCoincFilter5) static_cast<TH1*>(adummy->FindObject(Form("h_BS_PV_HF5_%s",varname[i])))->Fill(parameter);
          }
        }
        
        if (pprimaryVertexFilter) static_cast<TH1*>(adummy->FindObject(Form("h_PV_%s",varname[i])))->Fill(parameter);
        if (phfCoincFilter1) static_cast<TH1*>(adummy->FindObject(Form("h_HF1_%s",varname[i])))->Fill(parameter);
        if (phfCoincFilter2) static_cast<TH1*>(adummy->FindObject(Form("h_HF2_%s",varname[i])))->Fill(parameter);
        if (phfCoincFilter3) static_cast<TH1*>(adummy->FindObject(Form("h_HF3_%s",varname[i])))->Fill(parameter);
        if (phfCoincFilter4) static_cast<TH1*>(adummy->FindObject(Form("h_HF4_%s",varname[i])))->Fill(parameter);
        if (phfCoincFilter5) static_cast<TH1*>(adummy->FindObject(Form("h_HF5_%s",varname[i])))->Fill(parameter);
        
        
        if ((HLT_Trigger_part1 || HLT_Trigger_part2 || HLT_Trigger_part3 || HLT_Trigger_part4 || HLT_Trigger_part5 ||
             HLT_Trigger_part6 || HLT_Trigger_part7 || HLT_Trigger_part8 || HLT_Trigger_part9 || HLT_Trigger_part10 ||
             HLT_Trigger_part11 || HLT_Trigger_part12 || HLT_Trigger_part13 || HLT_Trigger_part14 || HLT_Trigger_part15 ||
             HLT_Trigger_part16 || HLT_Trigger_part17 || HLT_Trigger_part18 || HLT_Trigger_part19 || HLT_Trigger_part20))
        {
          
          static_cast<TH1*>(adummy->FindObject(Form("h_%s_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          static_cast<TH2*>(adummy->FindObject(Form("h_%s_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
          static_cast<TH2*>(adummy->FindObject(Form("h_%s_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
          static_cast<TH2*>(adummy->FindObject(Form("h_%s_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
          
          if (pBeamScrapingFilter)
          {
            static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
//            static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
//            static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
//            static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
            
            if (pprimaryVertexFilter)
            {
              static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
              static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
              static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
              static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
              
              if (phfCoincFilter1)
              {
                static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_HF1_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF1_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF1_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF1_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
              }
              if (phfCoincFilter2)
              {
                static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_HF2_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF2_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF2_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF2_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
              }
              if (phfCoincFilter3)
              {
                static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_HF3_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF3_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF3_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF3_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
              }
              else static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_NotHF3_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
              if (phfCoincFilter4)
              {
                static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_HF4_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF4_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF4_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF4_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
              }
              if (phfCoincFilter5)
              {
                static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_PV_HF5_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF5_hiHFasym_vs_%s",tLabel.Data(),varname[i])))->Fill(parameter,(hfplus-hfminus)/hf);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF5_hiHFasymLE_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingEAsym);
//                static_cast<TH2*>(adummy->FindObject(Form("h_%s_BS_PV_HF5_hiHFasymLET_vs_%s",tLabel.Data(),varN.Data())))->Fill(parameter,leadingETAsym);
              }
            }
            else
            {
              static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_NotPV_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
              
              if (phfCoincFilter3) static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_NotPV_HF3_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
              else static_cast<TH1*>(adummy->FindObject(Form("h_%s_BS_NotPV_NotHF3_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
            }
          }
          else static_cast<TH1*>(adummy->FindObject(Form("h_%s_NotBS_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          
          if (pprimaryVertexFilter) static_cast<TH1*>(adummy->FindObject(Form("h_%s_PV_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          else static_cast<TH1*>(adummy->FindObject(Form("h_%s_NotPV_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          
          if (phfCoincFilter1) static_cast<TH1*>(adummy->FindObject(Form("h_%s_HF1_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          else static_cast<TH1*>(adummy->FindObject(Form("h_%s_NotHF1_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          if (phfCoincFilter2) static_cast<TH1*>(adummy->FindObject(Form("h_%s_HF2_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          else static_cast<TH1*>(adummy->FindObject(Form("h_%s_NotHF2_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          if (phfCoincFilter3) static_cast<TH1*>(adummy->FindObject(Form("h_%s_HF3_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          else static_cast<TH1*>(adummy->FindObject(Form("h_%s_NotHF3_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          if (phfCoincFilter4) static_cast<TH1*>(adummy->FindObject(Form("h_%s_HF4_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
          else static_cast<TH1*>(adummy->FindObject(Form("h_%s_NotHF4_%s",tLabel.Data(),varN.Data())))->Fill(parameter);
        }
      }
    }
  }
  
  // save results
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/../../DataHistos";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  TFile *f = new TFile(Form("%s/distributions_%s_%s.root",saveDIR.c_str(),sampleLabel,tLabel.Data()),"RECREATE");
  aSave->Write("distributions", TObject::kOverwrite | TObject::kSingleKey);
  f->Close();
}

bool computeMCeff(const char* mcResults,
                  const char* hfFilter = "3",
                  const char* bsFilter = "BS",
                  const char* pvFilter = "PV",
                  const char* trigFilt = "MB")
{
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string histosDIR = mainDIR + "/../../DataHistos";
  
  TString effLabel(Form("%s_%s_%s_HF%s",trigFilt,bsFilter,pvFilter,hfFilter));
  if (!strcmp(trigFilt,"") && !strcmp(bsFilter,"") && !strcmp(hfFilter,"0") && !strcmp(pvFilter,"PV"))
  {
    effLabel.Clear();
    effLabel = "PV"; // FIXME: Other cases not implemented. Implement them!
  }
  cout << effLabel.Data() << endl;
  // Read file with mb distributions
  TFile *f = TFile::Open(Form("%s/distributions_%s_%s.root",histosDIR.c_str(),mcResults,!strcmp(trigFilt,"")?"MB":trigFilt));
  if (!f) {cout << "[ERROR] No file distributions_" << mcResults << "_" << trigFilt << ".root was found" << endl; return false;}
  
  // Get array with histos
  TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("distributions"));
  if (!arrHistos) {cout << "[ERROR] No histos array found in file " << f->GetName() << endl; return false;}
  //
  
  TObjArray* aSave = new TObjArray();
  
  TH1* h(0x0);
  TH1* h_MB_evtSel(0x0);
  
  const int nvar = 8;
  const char* varname[nvar] = {"hiHF", "hiHFplus","hiHFplusEta4", "hiHFminus","hiHFminusEta4","hiNtracks","hiHFhit","hiNpix"};
  
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    TObjArray* arrVar = static_cast<TObjArray*>(arrHistos->FindObject(Form("histos_%s",varN.Data())));
    if (!arrVar) {cout << "[ERROR] No histos array found for var " << varN.Data() << " found in file distributions_" << mcResults << ".root" << endl; return false;}
    
    h = static_cast<TH1*>(arrVar->FindObject(Form("h_%s",varN.Data()))->Clone());
    h_MB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("h_%s_%s",effLabel.Data(),varN.Data()))->Clone());
    
    if (!h || !h_MB_evtSel) {cout << "[ERROR] No histos found to compute efficiency " << endl; return false;}
    
    TEfficiency* eff = new TEfficiency(*h_MB_evtSel,*h);
    eff->SetName(Form("eff_%s_%s",effLabel.Data(),varN.Data()));
    aSave->Add(eff);
  }
  
  // save results
  TFile *fSave = new TFile(Form("%s/filterEff_%s_%s.root",histosDIR.c_str(),mcResults,effLabel.Data()),"RECREATE");
  aSave->Write("filterEff", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
  
  return true;
}

TEfficiency* getMCeff(const char* varname,
                      const char* mcResults,
                      const char* hfFilter = "3",
                      const char* bsFilter = "BS",
                      const char* pvFilter = "PV",
                      const char* trigFilt = "MB")
{
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string histosDIR = mainDIR + "/../../DataHistos";
  
  TString effLabel(Form("%s_%s_%s_HF%s",trigFilt,bsFilter,pvFilter,hfFilter));
  if (!strcmp(trigFilt,"") && !strcmp(bsFilter,"") && !strcmp(hfFilter,"0") && !strcmp(pvFilter,"PV"))
  {
    effLabel.Clear();
    effLabel = "PV"; // FIXME: Other cases not implemented. Implement them!
  }
  
  // Read file with eff distributions
  TFile *f = TFile::Open(Form("%s/filterEff_%s_%s.root",histosDIR.c_str(),mcResults,effLabel.Data()));
  if (!f) {
    cout << "[INFO] No file filterEff_" << mcResults << "_" << effLabel.Data() << ".root was found" << endl;
    cout << "[INFO] Computing efficiencies for " << effLabel.Data() << " selection" << endl;
    if (!computeMCeff(mcResults,hfFilter,bsFilter,pvFilter,trigFilt)) {cout << "[ERROR] Efficiency could not be computed" << endl; return 0x0;}
    f = TFile::Open(Form("%s/filterEff_%s_%s.root",histosDIR.c_str(),mcResults,effLabel.Data()));
  }
  
  TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("filterEff"));
  if (!arrHistos) {cout << "[ERROR] No efficiency array found in file " << f->GetName() << endl; return 0x0;}
  
  TEfficiency* eff = static_cast<TEfficiency*>(arrHistos->FindObject(Form("eff_%s_%s",effLabel.Data(),varname))->Clone());
  f->Close();
  if (!eff) {cout << "[ERROR] No efficiency " << Form("eff_%s_%s",effLabel.Data(),varname) << " was found" << endl; return 0x0;}
  else {return eff;}
}
