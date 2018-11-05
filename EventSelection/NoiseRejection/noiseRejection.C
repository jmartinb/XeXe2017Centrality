// computes level of noise rejection by different event filters
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
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
#include <TCut.h>
//canvas, legend, latex ... //cosmetic
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>
//random
#include <TRandom.h>
#include <TStopwatch.h>
#include <ctime>
//private setup
#include "../../Utilities/drawUtils.h"
#include "../mbDistributions.C"

void noiseRejection(bool recreateDistributions = false,
                    const char* sampleLabel = "DATA_EB",
                    const char* inputFile = "ebInputFiles.txt",
                    const char* triggerName = "HLT_HIL1NotBptxOR_v1",
                    const char* colLabel = "XeXe 5.44 TeV"
                    )
{
  
  // Run macro to create the noise distributions
  if (recreateDistributions) mbDistributions(Form("%s/../InputFiles/%s",gSystem->pwd(),inputFile),triggerName,sampleLabel);
  //
  
  TString tLabel(triggerName);
  if (tLabel.Contains("UnpairedBunches")) tLabel = "UB";
  else if (tLabel.Contains("NotBptxOR")) tLabel = "NBOR";
  else {cout << "[ERROR] Unkown trigger label" << endl; return;}
  
  // Read file with noise distributions
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string resultsDIR = mainDIR + "/Results/";
  string histosDIR = mainDIR + "/../../DataHistos";
  TFile *f = TFile::Open(Form("%s/distributions_%s_%s.root",histosDIR.c_str(),sampleLabel,tLabel.Data()));
  if (!f) {cout << "[ERROR] No file distributions_" << sampleLabel << "_" << tLabel.Data() << ".root was found" << endl; return;}
  //
  
  // Get array with histos
  TObjArray* arrHistosEB = static_cast<TObjArray*>(f->FindObjectAny("distributions"));
  if (!arrHistosEB) {cout << "ERROR: No histos array found in file " << f << endl; return;}
  //
  
  TObjArray* aSave = new TObjArray();
  
  const int nvar = 8;
  const char* varname[8] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus", "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix"};
  
  // Create directory to save results
  string saveDIR = mainDIR + "/Results/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  // Make the plots
  TObjArray* arrVarEB(0x0);
  TH1* hMB(0x0);
  TH1* hMB_PV(0x0);
  TH1* hMB_BS(0x0);
  TH1* hMB_CC(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    
    arrVarEB = static_cast<TObjArray*>(arrHistosEB->FindObject(Form("histos_%s",varN.Data())));
    
    hMB = static_cast<TH1*>(arrVarEB->FindObject(Form("h_%s_%s",tLabel.Data(),varN.Data()))->Clone());
    hMB_PV = static_cast<TH1*>(arrVarEB->FindObject(Form("h_%s_PV_%s",tLabel.Data(),varN.Data()))->Clone());
    hMB_BS = static_cast<TH1*>(arrVarEB->FindObject(Form("h_%s_BS_%s",tLabel.Data(),varN.Data()))->Clone());
    hMB_CC = static_cast<TH1*>(arrVarEB->FindObject(Form("h_%s_HF3_%s",tLabel.Data(),varN.Data()))->Clone());
    
    if (!hMB || !hMB_PV || !hMB_BS || !hMB_CC) {cout << "[ERROR]: No histos found in array for variable " << varN.Data() << endl; return;}
    
    double nEB = hMB->GetEntries();
    double nEB_PV = hMB_PV->GetEntries();
    double nEB_BS = hMB_BS->GetEntries();
    double nEB_CC = hMB_CC->GetEntries();
    
    double rej_PV = ((nEB-nEB_PV)/nEB)*100.;
    double rej_BS = ((nEB-nEB_BS)/nEB)*100.;
    double rej_CC = ((nEB-nEB_CC)/nEB)*100.;
    
    gStyle -> SetOptStat(0);
    TCanvas* c=  new TCanvas(Form("c_%s",varN.Data()),"", 900,500);
    c->Divide(2,1);
    c->cd(1);
    gPad->SetLogy();
    hMB->GetXaxis()->SetRangeUser(0,20);
    hMB->GetYaxis()->SetRangeUser(1,hMB->GetMaximum()*100.);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      hMB->SetXTitle (Form("%s #sum{E_T} (GeV)",axis.Data()));
    }
    hMB->SetLineColor(1);
    hMB->SetMarkerColor(1);
    hMB->SetMarkerStyle(20);
    hMB->DrawCopy();

    TLegend* l1 = new TLegend(0.098,0.77,0.7,0.80);
    legStyle(l1);
    
    l1->AddEntry(hMB, triggerName, "p");
    l1->Draw("same");
    drawText(colLabel,0.2,0.2+0.65);
    
    c->cd(2);
    if (hMB_PV->GetMaximum()>0) gPad->SetLogy();
    hMB_PV->GetXaxis()->SetRangeUser(0,20);
    hMB_PV->GetYaxis()->SetRangeUser(0,(hMB_PV->GetMaximum()>0)?(hMB_PV->GetMaximum()*100.):1.);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      hMB_PV->SetXTitle (Form("%s #sum{E_T} (GeV)",axis.Data()));
    }
    hMB_PV->SetLineColor(1);
    hMB_PV->SetMarkerColor(1);
    hMB_PV->SetMarkerStyle(20);
    hMB_PV->Draw();
    TLegend* l2 = new TLegend(0.098,0.7,0.7,0.73);
    legStyle(l2);
    l2->AddEntry(hMB_PV, Form("%s + vtx filter",triggerName), "p");
    drawText(Form("vtx filter rejected %s = %0.3f %%",triggerName,rej_PV),0.15,0.3);
    drawText(Form("bs filter rejected %s = %0.3f %%",triggerName,rej_BS),0.15,0.25);
    drawText(Form("HF3 filter rejected %s = %0.3f %%",triggerName,rej_CC),0.15,0.2);
    
    l2->Draw("same");
    
    c->SaveAs(Form("%s/noiseRejection_%s.pdf",saveDIR.c_str(),varN.Data()));
    aSave->Add(c);
  }
  
  // save results
  TFile *fSave = new TFile(Form("%s/noiseRejection_%s.root",saveDIR.c_str(),tLabel.Data()),"RECREATE");
  aSave->Write("noiseRejection", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
  f->Close();
  
  return;
}
