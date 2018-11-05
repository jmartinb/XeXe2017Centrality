// compute fraction of EM events after event selection
// Author : Javier Martin Blanco 25/01/2017

//basic c++ header, string ...
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <cstring>
#include <vector>
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
#include <TF1.h>
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

//double L = 3423000; //b^-1
double truncEt1_high = 36.8;
double truncEt1_low = 0.;
const char* centPerc1 = "80-100"; //Centrality percentile corresponding to above HF cut
double truncEt2_high = 78.6;
double truncEt2_low = 17.1;
const char* centPerc2 = "70-90"; //Centrality percentile corresponding to above HF cut
double truncEt3_high = 36.8;
double truncEt3_low = 24.8;
const char* centPerc3 = "80-85"; //Centrality percentile corresponding to above HF cut
double truncEt4_high = 24.8;
double truncEt4_low = 17.1;
const char* centPerc4 = "85-90"; //Centrality percentile corresponding to above HF cut

void emContamination(const char* hfFilter = "3", // Number of required towers in HF coincidence filter
                     const std::vector<bool> recreateDistributions = {false,false,false}, // {data,MC,EM}
                     const std::vector<const char*> inputFile = {"dataMBInputFiles.txt","emInputFiles_single.txt","emInputFiles_double.txt"},
                     const std::vector<const char*> sampleLabel = {"DATA_MB","STARlight_single","STARlight_double"},
                     const std::vector<double> xSection = {5.61,10.0,0.041}, // cross section: {sigma_had,sigma_EM_single,sigma_EM_double} b
                     const std::vector<double> eff = {0.95,0.0044,0.3226}, // selection eff: {eff+cont_data,eff_EM_single,eff_EM_double}
                     const char* colLabel = "XeXe 5.44 TeV",
                     const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1"
                 )
{
  
  TString tLabel(triggerName);
  if (tLabel.Contains("MinimumBias")) tLabel = "MB";
  else {cout << "[ERROR] Unkown trigger label" << endl; return;}
  
  // Check the input parameters
  const int nSamples = 3;
  if (strcmp(hfFilter,"1") && strcmp(hfFilter,"2") && strcmp(hfFilter,"3") && strcmp(hfFilter,"4")) {cout << "[ERROR] hfFilter must be the number of towers for the filter. Options are: 1 , 2 , 3 or 4 " << endl; return;}
  if (recreateDistributions.size()!=nSamples) {cout << "[ERROR] recreateDistributions must have " << nSamples << " entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  if (inputFile.size()!=nSamples) {cout << "[ERROR] inputfile must have " << nSamples << " entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  if (sampleLabel.size()!=nSamples) {cout << "[ERROR] sampleLabel must have " << nSamples << " entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  if (xSection.size()!=nSamples) {cout << "[ERROR] xSection must have 3 entries. Usage: {data,starlight_single,starlight_double}" << endl; return;}
  
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string resultsDIR = mainDIR + "/Results";
  string histosDIR = mainDIR + "/../../DataHistos";
  
  TObjArray* arrHistosarr = new TObjArray();
  for (int n = 0 ; n < nSamples ; n++)
  {
    // Run macro to create the mb distributions
    if (recreateDistributions[n])
    {
      cout << "[INFO] Creating distributions for " << sampleLabel[n]  << endl;
      mbDistributions(Form("%s/../../InputFiles/%s",gSystem->pwd(),inputFile[n]),triggerName,sampleLabel[n]);
    }
    
    // Read file with mb distributions
    TFile *f = TFile::Open(Form("%s/distributions_%s_%s.root",histosDIR.c_str(),sampleLabel[n],tLabel.Data()));
    if (!f) {cout << "[ERROR] No file distributions_" << sampleLabel[n] << "_" << tLabel.Data() << ".root was found" << endl; return;}
    //
    // Get array with histos
    TObjArray* arrHistos = static_cast<TObjArray*>(f->FindObjectAny("distributions"));
    if (!arrHistos) {cout << "ERROR: No histos array found in file " << f->GetName() << endl; return;}
    arrHistosarr->Add(arrHistos);
    //
  }
  
  TObjArray* aSave = new TObjArray();
  aSave->SetOwner(true);
  
  const int nvar = 10;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus",
    "hiHFminusEta4","hiNtracks","hiHFhit","hiNpix","nTrk","nTrkVtx"};
  
  // Create directory to save the plots
  string saveDIR = mainDIR + "/Results/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  // Define cross section weigthing functions
  TF1* fW = new TF1("wEMS","[0]",0.,1000.);
  
  // Make the plots
  TObjArray* arrVar(0x0);
  TH1* hMB_evtSel(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    
    gStyle -> SetOptStat(0);
    TCanvas* c=  new TCanvas(Form("c_%s",varN.Data()),"", 900,500);
    c->Divide(2,1);
    
    TLegend* l1 = new TLegend(0.45,0.64,1.05,0.76);
    legStyle(l1);
    
    double nMB_evtSel(0.),nEMS_evtSel(0.),nEMD_evtSel(0.);
    double nMB_evtSel_trunc1(0.),nEMS_evtSel_trunc1(0.),nEMD_evtSel_trunc1(0.);
    double nMB_evtSel_trunc2(0.),nEMS_evtSel_trunc2(0.),nEMD_evtSel_trunc2(0.);
    double nMB_evtSel_trunc3(0.),nEMS_evtSel_trunc3(0.),nEMD_evtSel_trunc3(0.);
    double nMB_evtSel_trunc4(0.),nEMS_evtSel_trunc4(0.),nEMD_evtSel_trunc4(0.);
    double wEMS(0.),wEMD(0.);//,wMB(0.);
//    double sigma_mb = xSection[0] + xSection[1] + xSection[2];
    
    TObjArray* arrHistosSave = new TObjArray();
    arrHistosSave->SetOwner(true);
    
    for (int n = 0 ; n < nSamples ; n++)
    {
      TObjArray* arrHistos = static_cast<TObjArray*>(arrHistosarr->At(n));
      if (arrHistos) arrVar = static_cast<TObjArray*>(arrHistos->FindObject(Form("histos_%s",varN.Data())));
      else {cout << "[ERROR] No array of histos found for " << sampleLabel[n] << endl; return;}
      
      hMB_evtSel = static_cast<TH1*>(arrVar->FindObject(Form("h_MB_BS_PV_HF%s_%s",hfFilter,varN.Data()))->Clone(Form("h_MB_BS_PV_HF%s_%s_EMsubs",hfFilter,varN.Data())));
      hMB_evtSel->Sumw2();
      if (!hMB_evtSel)
      {
        cout << "[ERROR] No MB histos found in array " << endl;
        return;
      }
      
      if (n==0)
      {
        nMB_evtSel = hMB_evtSel->GetEntries();
        nMB_evtSel_trunc1 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt1_low),hMB_evtSel->FindBin(truncEt1_high));
        nMB_evtSel_trunc2 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt2_low),hMB_evtSel->FindBin(truncEt2_high));
        nMB_evtSel_trunc3 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt3_low),hMB_evtSel->FindBin(truncEt3_high));
        nMB_evtSel_trunc4 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt4_low),hMB_evtSel->FindBin(truncEt4_high));
//        wMB = (L*sigma_mb*eff[0])/nMB;
//        fW->SetParameter(0,wMB);
//        hMB_evtSel->Multiply(fW);
      }
      else if (n==1)
      {
        nEMS_evtSel = hMB_evtSel->GetEntries();
        nEMS_evtSel_trunc1 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt1_low),hMB_evtSel->FindBin(truncEt1_high));
        nEMS_evtSel_trunc2 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt2_low),hMB_evtSel->FindBin(truncEt2_high));
        nEMS_evtSel_trunc3 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt3_low),hMB_evtSel->FindBin(truncEt3_high));
        nEMS_evtSel_trunc4 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt4_low),hMB_evtSel->FindBin(truncEt4_high));
        wEMS = (nMB_evtSel/nEMS_evtSel)*(eff[1]*xSection[1])/(eff[0]*xSection[0] + eff[1]*xSection[1] + eff[2]*xSection[2]);//(nMB/sigma_mb)*(xSection[1]/nEMS);
        fW->SetParameter(0,wEMS);
        hMB_evtSel->Multiply(fW);
      }
      else if (n==2)
      {
        nEMD_evtSel = hMB_evtSel->GetEntries();
        nEMD_evtSel_trunc1 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt1_low),hMB_evtSel->FindBin(truncEt1_high));
        nEMD_evtSel_trunc2 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt2_low),hMB_evtSel->FindBin(truncEt2_high));
        nEMD_evtSel_trunc3 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt3_low),hMB_evtSel->FindBin(truncEt3_high));
        nEMD_evtSel_trunc4 = hMB_evtSel->Integral(hMB_evtSel->FindBin(truncEt4_low),hMB_evtSel->FindBin(truncEt4_high));
        wEMD = (nMB_evtSel/nEMD_evtSel)*(eff[2]*xSection[2])/(eff[0]*xSection[0] + eff[1]*xSection[1] + eff[2]*xSection[2]);//(nMB/sigma_mb)*(xSection[2]/nEMD);
        fW->SetParameter(0,wEMD);
        hMB_evtSel->Multiply(fW);
      }
      else {cout << "[ERROR] There are more samples than expected" << endl; return;}
      
      c->cd(1);
//      gPad->SetLogy();
      hMB_evtSel->GetXaxis()->SetRangeUser(0,50);
      hMB_evtSel->SetLineColor(n+1);
      hMB_evtSel->SetMarkerColor(n+1);
      hMB_evtSel->SetMarkerStyle(20);
      if(varN.Contains("hiHF") && !varN.Contains("hit"))
      {
        TString axis(varN.Data());
        axis.Replace(0,2,"");
        hMB_evtSel->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
        hMB_evtSel->GetXaxis()->SetTitleOffset(1.3);
      }
      hMB_evtSel->DrawClone(n==0?"":"same");
      
      l1->AddEntry(hMB_evtSel, sampleLabel[n], "p");
      
      arrHistosSave->Add(hMB_evtSel);
    }
    
    
    // Compute EM contamination
    double emCont = (eff[1]*xSection[1] + eff[2]*xSection[2])/(eff[0]*xSection[0] + eff[1]*xSection[1] + eff[2]*xSection[2])*100.;
    double emCont_subs = ((wEMS*nEMS_evtSel)+(wEMD*nEMD_evtSel))*100./(nMB_evtSel);
    double emCont_subs_ETtrunc1 = ((wEMS*nEMS_evtSel_trunc1)+(wEMD*nEMD_evtSel_trunc1))*100./(nMB_evtSel_trunc1);
    double emCont_subs_ETtrunc2 = ((wEMS*nEMS_evtSel_trunc2)+(wEMD*nEMD_evtSel_trunc2))*100./(nMB_evtSel_trunc2);
    double emCont_subs_ETtrunc3 = ((wEMS*nEMS_evtSel_trunc3)+(wEMD*nEMD_evtSel_trunc3))*100./(nMB_evtSel_trunc3);
    double emCont_subs_ETtrunc4 = ((wEMS*nEMS_evtSel_trunc4)+(wEMD*nEMD_evtSel_trunc4))*100./(nMB_evtSel_trunc4);
    
    // EM contamination subtraction from MB
    TH1* hMB_subs = static_cast<TH1*>(arrHistosSave->At(0)->Clone());
    TH1* hEMs = static_cast<TH1*>(arrHistosSave->At(1));
    TH1* hEMd = static_cast<TH1*>(arrHistosSave->At(2));
    hMB_subs->Add(hEMs,-1.);
    hMB_subs->Add(hEMd,-1.);
    hMB_subs->SetLineColor(4);
    hMB_subs->SetMarkerColor(4);
    hMB_subs->SetMarkerStyle(21);
    l1->AddEntry(hMB_subs, "DATA (EM-subtracted)", "p");
    
    c->cd(1);
    hMB_subs->DrawClone("same");
    l1->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(Form("MB + vtx +BS + HF%s filters",hfFilter),0.54,0.2+0.60);

    // EM contamination vs var
    TH1* hMB = static_cast<TH1*>(arrHistosSave->At(0));
    TH1* hEM_tot = static_cast<TH1*>(hEMs->Clone("EM_tot"));
    hEM_tot->Add(hEMd);
    hEM_tot->Divide(hMB);
    hEM_tot->GetXaxis()->SetRangeUser(0,50);
    hEM_tot->SetYTitle("% EM cont. /100");
    hEM_tot->GetYaxis()->SetTitleOffset(1.6);
    hEM_tot->GetYaxis()->SetRangeUser(0.,0.35);
    hEM_tot->SetLineColor(1);
    hEM_tot->SetMarkerColor(1);
    hEM_tot->SetMarkerStyle(22);
    
    c->cd(2);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      hEM_tot->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
    }
    hEM_tot->Draw();
    drawText(Form("EM contamination vs %s",varN.Data()),0.54,0.2+0.65);
    drawText(Form("Global EM cont. = %0.2f %%",emCont_subs),0.40,0.80);
    if(!varN.CompareTo("hiHF"))
    {
      drawText(Form("EM cont. %s%% (%0.1f < #sum E_{T} < %0.1f) = %0.2f %%",centPerc1,truncEt1_low,truncEt1_high,emCont_subs_ETtrunc1),0.18,0.75);
      drawText(Form("EM cont. %s%% (%0.1f < #sum E_{T} < %0.1f) = %0.2f %%",centPerc2,truncEt2_low,truncEt2_high,emCont_subs_ETtrunc2),0.18,0.70);
      drawText(Form("EM cont. %s%% (%0.1f < #sum E_{T} < %0.1f) = %0.2f %%",centPerc3,truncEt3_low,truncEt3_high,emCont_subs_ETtrunc3),0.18,0.65);
      drawText(Form("EM cont. %s%% (%0.1f < #sum E_{T} < %0.1f) = %0.2f %%",centPerc4,truncEt4_low,truncEt4_high,emCont_subs_ETtrunc4),0.18,0.60);
    }
    
    c->SaveAs(Form("%s/emContamination_%s_HF%s.pdf",saveDIR.c_str(),varN.Data(),hfFilter));
    aSave->Add(c);
    aSave->Add(hMB_subs);
    
//    delete arrHistosSave;
  }

  // save results
  TFile *fSave = new TFile(Form("%s/emContamination_HF%s.root",saveDIR.c_str(),hfFilter),"RECREATE");
  aSave->Write("emContamination", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
}
