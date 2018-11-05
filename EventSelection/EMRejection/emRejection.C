// compute fraction of EM events rejected by event selection
// Author : Javier Martin Blanco 25/01/2017

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

void emRejection(bool recreateDistributions = false,
                     const char* inputFile = "emInputFiles_single.txt",
                     const char* sampleLabel = "STARlight_single",
                     const char* colLabel = "XeXe 5.44 TeV",
                     const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1"
                 )
{
  TString tLabel(triggerName);
  if (tLabel.Contains("MinimumBias")) tLabel = "MB";
  else {cout << "[ERROR] Unkown trigger label" << endl; return;}
  
  // Run macro to create the mb distributions
  if (recreateDistributions) mbDistributions(Form("%s/../../InputFiles/%s",gSystem->pwd(),inputFile),triggerName,sampleLabel);
  //
  
  // Read file with mb distributions
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string histosDIR = mainDIR + "/../../DataHistos";
  string resultsDIR = mainDIR + "/Results";
  TFile *f = TFile::Open(Form("%s/distributions_%s_%s.root",histosDIR.c_str(),sampleLabel,tLabel.Data()));
  if (!f) {cout << "[ERROR] No file distributions_" << sampleLabel << "_" << tLabel.Data() << ".root was found" << endl; return;}
  //
  
  // Get array with histos
  TObjArray* arrHistosMB = static_cast<TObjArray*>(f->FindObjectAny("distributions"));
  if (!arrHistosMB) {cout << "[ERROR] No histos array found in file " << f << endl; return;}
  //
  
  TObjArray* aSave = new TObjArray();
  
  const int nvar = 12;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus",
    "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix","nTrk","trkEta","trkPt","nTrkVtx"};
  
  // Create directory to save the plots
  TString label(sampleLabel);
  string saveDIR = mainDIR + "/Results";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  // Make the plots
  TObjArray* arrVar(0x0);
  TH1* h(0x0);
  TH1* hMB(0x0);
  TH1* hMB_noNoise(0x0);
  TH1* hMB_evtSel(0x0);
  TH1* hMB_evtSel2(0x0);
  TH1* hMB_evtSel3(0x0);
  TH1* hMB_evtSel4(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    
    arrVar = static_cast<TObjArray*>(arrHistosMB->FindObject(Form("histos_%s",varN.Data())));
    
    h = static_cast<TH1*>(arrVar->FindObject(Form("h_%s",varN.Data()))->Clone());
    hMB = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_%s",varN.Data()))->Clone());
    hMB_noNoise = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_BS_PV_%s",varN.Data()))->Clone());
    hMB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_BS_PV_HF1_%s",varN.Data()))->Clone());
    hMB_evtSel2 = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_BS_PV_HF2_%s",varN.Data()))->Clone());
    hMB_evtSel3 = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_BS_PV_HF3_%s",varN.Data()))->Clone());
    hMB_evtSel4 = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_BS_PV_HF4_%s",varN.Data()))->Clone());
    if (!h || !hMB || !hMB_noNoise || !hMB_evtSel || !hMB_evtSel2 || !hMB_evtSel3 || !hMB_evtSel4)
    {
      cout << "[ERROR] No MB histos found in array " << endl;
      return;
    }
    
    double n = h->GetEntries();
    double nMB = hMB->GetEntries();
    double nMB_noNoise = hMB_noNoise->GetEntries();
    double nMB_evtSel = hMB_evtSel->GetEntries();
    double nMB_evtSel2 = hMB_evtSel2->GetEntries();
    double nMB_evtSel3 = hMB_evtSel3->GetEntries();
    double nMB_evtSel4 = hMB_evtSel4->GetEntries();
    
    double rej_trigger = ((n-nMB)/n)*100.;
    double rej_noNoise = ((n-nMB_noNoise)/n)*100.;
    double rej_evtSel = ((n-nMB_evtSel)/n)*100.;
    double rej_evtSel2 = ((n-nMB_evtSel2)/n)*100.;
    double rej_evtSel3 = ((n-nMB_evtSel3)/n)*100.;
    double rej_evtSel4 = ((n-nMB_evtSel4)/n)*100.;
    
    double xMin = 0.;
    double xMax = 50.;
    bool logY = false;
    if(!strcmp(varN.Data(),"trkEta") ) {xMin = -4.; xMax = 4.; logY = true;}
    if(!strcmp(varN.Data(),"trkPt") ) {xMin = 0.; xMax = 15.; logY = true;}
    
    gStyle -> SetOptStat(0);
    TCanvas* c=  new TCanvas(Form("c_%s",varN.Data()),"", 900,500);
    c->Divide(2,1);
    c->cd(1);
    if (logY) gPad->SetLogy();
    h->GetXaxis()->SetRangeUser(xMin,xMax);
    h->SetLineColor(1);
    h->SetMarkerColor(1);
    h->SetMarkerStyle(21);
    h->DrawCopy();
    
    hMB->SetLineColor(2);
    hMB->SetMarkerColor(2);
    hMB->SetMarkerStyle(20);
    hMB->DrawCopy("same");
    
    TLegend* l1 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l1);
    
    l1->AddEntry(h, "All EM events", "p");
    l1->AddEntry(hMB, "MB selected EM events", "p");
    l1->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(sampleLabel,0.54,0.2+0.60);
    drawText(Form("MB trigger rejected = %0.3f %%",rej_trigger),0.15,0.35);
    
    c->cd(2);
    if (logY) gPad->SetLogy();
    hMB_noNoise->GetXaxis()->SetRangeUser(xMin,xMax);
    hMB_noNoise->SetLineColor(3);
    hMB_noNoise->SetMarkerColor(3);
    hMB_noNoise->SetMarkerStyle(20);
    hMB_noNoise->Draw();
    //
    hMB_evtSel->GetXaxis()->SetRangeUser(xMin,xMax);
    hMB_evtSel->SetLineColor(4);
    hMB_evtSel->SetMarkerColor(4);
    hMB_evtSel->SetMarkerStyle(20);
    hMB_evtSel->Draw("same");
    //
    hMB_evtSel2->GetXaxis()->SetRangeUser(xMin,xMax);
    hMB_evtSel2->SetLineColor(5);
    hMB_evtSel2->SetMarkerColor(5);
    hMB_evtSel2->SetMarkerStyle(20);
    hMB_evtSel2->Draw("same");
    //
    hMB_evtSel3->GetXaxis()->SetRangeUser(xMin,xMax);
    hMB_evtSel3->SetLineColor(6);
    hMB_evtSel3->SetMarkerColor(6);
    hMB_evtSel3->SetMarkerStyle(20);
    hMB_evtSel3->Draw("same");
    //
    hMB_evtSel4->GetXaxis()->SetRangeUser(xMin,xMax);
    hMB_evtSel4->SetLineColor(7);
    hMB_evtSel4->SetMarkerColor(7);
    hMB_evtSel4->SetMarkerStyle(20);
    hMB_evtSel4->Draw("same");
    //
    
    TLegend* l2 = new TLegend(0.37,0.55,0.97,0.81);
    legStyle(l2);
    l2->AddEntry(hMB_noNoise,"MB + vtx +BS filters", "p");
    l2->AddEntry(hMB_evtSel, "MB + vtx +BS + HF1 filters", "p");
    l2->AddEntry(hMB_evtSel2, "MB + vtx +BS + HF2 filters", "p");
    l2->AddEntry(hMB_evtSel3, "MB + vtx +BS + HF3 filters", "p");
    l2->AddEntry(hMB_evtSel4, "MB + vtx +BS + HF4 filters", "p");
    drawText(Form("MB + vtx + BS filters rejected = %0.3f %%",rej_noNoise),0.15,0.35);
    drawText(Form("MB + vtx + BS + HF1 filters rejected = %0.3f %%",rej_evtSel),0.15,0.3);
    drawText(Form("MB + vtx + BS + HF2 filters rejected = %0.3f %%",rej_evtSel2),0.15,0.25);
    drawText(Form("MB + vtx + BS + HF3 filters rejected = %0.3f %%",rej_evtSel3),0.15,0.20);
    drawText(Form("MB + vtx + BS + HF4 filters rejected = %0.3f %%",rej_evtSel4),0.15,0.15);
    
    l2->Draw("same");
    
    c->SaveAs(Form("%s/emRejection_%s_%s.pdf",saveDIR.c_str(),varN.Data(),label.Contains("single")?"single":(label.Contains("double")?"double":label.Data())));
    aSave->Add(c);
  }

  // save results
  TFile *fSave = new TFile(Form("%s/emRejection_%s.root",saveDIR.c_str(),label.Contains("single")?"single":(label.Contains("double")?"double":label.Data())),"RECREATE");
  aSave->Write("emRejection", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
  f->Close();
}
