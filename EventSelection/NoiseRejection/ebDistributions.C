// Produce empty bunches distributions with different event filters applied
// Author : Javier Martin Blanco 24/01/2018

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH1D.h>
#include <TObjArray.h>

#include "../../Utilities/treeUtils.h"


//////////////////////////////////////////////
//////////////////////////////////////////////

void ebDistributions(bool recreateDistributions = false,
                     const char* sampleLabel = "DATA_EB",
                     const char* inputFile = "ebInputFiles.txt",
                     const char* triggerName = "HLT_HINoBptxOR_v1"
                     )
{
  if (recreateDistributions)
  {
    cout << "[INFO] Creating distributions for " << sampleLabel  << endl;
    mbDistributions(inputFile,triggerName,sampleLabel);
  }
  
  // Read files and chain TTrees
  TChain* tree = readFiles(inputFile);
  //
  
  unsigned int nEvents = tree->GetEntries();
  if (nEvents>0) cout << "[INFO] Number of events to be processed = " << nEvents << endl;
  else {cout << "[ERROR] No events to be processed" << endl; return;}
  
  // Load needed branches
  if (!loadBranches(triggerName,tree)) {cout << "[ERROR] Error loading branches" << endl; return;}
  //

  // Definition of variable ranges and number of bins for plots
  double hiHFMax = 300;
  double hiHFSSMax = 200;
  double hiHFSSTruncMax = 80;
  double hiHFhitMax = 8000;
  double hiNpixMax = 2000;
  double hiNtracksMax = 400;
  
  int nBin = 100;
  //
  
  // Definition of histos to store distributions
  const int nvar = 8;
  const char* varname[8] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus", "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix"};
  
  TObjArray* aSave = new TObjArray();
  
  for (int i = 0 ; i < nvar ; i++)
  {
    TObjArray* arr = new TObjArray();
    arr->SetName(Form("histos_%s",varname[i]));
    
    double xMin = 0.;
    double xMax = 300;
    
    if( !strcmp(varname[i],"hiNtracks") ) {xMax = hiNtracksMax; nBin = (int)hiNtracksMax;}
    if( !strcmp(varname[i],"hiHF")) {xMax = hiHFMax; nBin = (int)hiHFMax*8;}
    if( !strcmp(varname[i],"hiHFplus") || !strcmp(varname[i],"hiHFminus")) {xMax = hiHFSSMax; nBin = (int)hiHFSSMax*8;}
    if( !strcmp(varname[i],"hiHFplusEta4") || !strcmp(varname[i],"hiHFminusEta4")) {xMax = hiHFSSTruncMax; nBin = (int)hiHFSSTruncMax*8;}
    if( !strcmp(varname[i],"hiHFhit") ){xMax = hiHFhitMax; nBin = hiHFhitMax;}
    if( !strcmp(varname[i],"hiNpix") ){xMax = hiNpixMax; nBin = hiNpixMax;}
    
    TH1D* hEB = new TH1D(Form("hEB_%s",varname[i]), Form(";%s;",varname[i]), nBin,xMin,xMax);
    TH1D* hEB_BS = (TH1D*)hEB->Clone(Form("hEB_BS_%s",varname[i]));
    TH1D* hEB_CC = (TH1D*)hEB->Clone(Form("hEB_CC_%s",varname[i]));
    TH1D* hEB_PV = (TH1D*)hEB->Clone(Form("hEB_PV_%s",varname[i]));
    TH1D* hEB_HF = (TH1D*)hEB->Clone(Form("hEB_HF_%s",varname[i]));
    TH1D* hEB_BS_PV = (TH1D*)hEB->Clone(Form("hEB_BS_PV_%s",varname[i]));
    TH1D* hEB_BS_PV_HF1 = (TH1D*)hEB->Clone(Form("hEB_BS_PV_HF1_%s",varname[i]));
    TH1D* hEB_CC_PV = (TH1D*)hEB->Clone(Form("hEB_CC_PV_%s",varname[i]));
    TH1D* hEB_CC_PV_HF1 = (TH1D*)hEB->Clone(Form("hEB_CC_PV_HF1_%s",varname[i]));
    
    arr->Add(hEB);
    arr->Add(hEB_BS);
    arr->Add(hEB_CC);
    arr->Add(hEB_PV);
    arr->Add(hEB_HF);
    arr->Add(hEB_BS_PV);
    arr->Add(hEB_BS_PV_HF1);
    arr->Add(hEB_CC_PV);
    arr->Add(hEB_CC_PV_HF1);
    
    aSave->Add(arr);
  }
  //
  
  TObjArray* adummy(0x0);
  // Event loop
  for(unsigned int iev = 0; iev < nEvents; iev++)
  {
    if(iev%5000 == 0) cout<<"Processing event: " << iev << " / " << nEvents << endl;
    tree->GetEntry(iev);
  
    for (int i = 0 ; i < nvar ; i++)
    {
      adummy = static_cast<TObjArray*>(aSave->FindObject(Form("histos_%s",varname[i])));
      
      float parameter = -1;
      if(!strcmp(varname[i],"hiHF")) parameter = hf;
      if(!strcmp(varname[i],"hiHFplus")) parameter = hfplus;
      if(!strcmp(varname[i],"hiHFminus")) parameter = hfminus;
      if(!strcmp(varname[i],"hiHFplusEta4")) parameter = hfpluseta4;
      if(!strcmp(varname[i],"hiHFminusEta4")) parameter = hfminuseta4;
      if(!strcmp(varname[i],"hiNtracks")) parameter = hiNtrks;
      if(!strcmp(varname[i],"hiHFhit")) parameter = hfhit;
      if(!strcmp(varname[i],"hiNpix")) parameter = npix;
      
      if (HLT_Trigger_part1)
      {
        static_cast<TH1*>(adummy->FindObject(Form("hEB_%s",varname[i])))->Fill(parameter);
        
        if (pBeamScrapingFilter)
        {
          static_cast<TH1*>(adummy->FindObject(Form("hEB_BS_%s",varname[i])))->Fill(parameter);
          
          if (pprimaryVertexFilter)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hEB_BS_PV_%s",varname[i])))->Fill(parameter);
            
            if (phfCoincFilter1)
            {
              static_cast<TH1*>(adummy->FindObject(Form("hEB_BS_PV_HF1_%s",varname[i])))->Fill(parameter);
            }
          }
        }
        
        if (pclusterCompatibilityFilter)
        {
          static_cast<TH1*>(adummy->FindObject(Form("hEB_CC_%s",varname[i])))->Fill(parameter);
          
          if (pprimaryVertexFilter)
          {
            static_cast<TH1*>(adummy->FindObject(Form("hEB_CC_PV_%s",varname[i])))->Fill(parameter);
            
            if (phfCoincFilter1)
            {
              static_cast<TH1*>(adummy->FindObject(Form("hEB_CC_PV_HF1_%s",varname[i])))->Fill(parameter);
            }
          }
        }
        
        if (pprimaryVertexFilter)
        {
          static_cast<TH1*>(adummy->FindObject(Form("hEB_PV_%s",varname[i])))->Fill(parameter);
        }
        
        if (phfCoincFilter1)
        {
          static_cast<TH1*>(adummy->FindObject(Form("hEB_HF_%s",varname[i])))->Fill(parameter);
        }
      }
    }
  }
  
  // save results
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/Results/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  TFile *f = new TFile(Form("%s/ebDistributions.root",saveDIR.c_str()),"RECREATE");
  aSave->Write("ebDistributions", TObject::kOverwrite | TObject::kSingleKey);
  f->Close();
  
  return;
}
