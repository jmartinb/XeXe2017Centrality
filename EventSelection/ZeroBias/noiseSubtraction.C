// Subtraction of noise from ZB sample
// This macro uses the normalisation wrt to the total number of BX in the data taking
// Author : Javier Martin Blanco 04/04/2018

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
#include <TEfficiency.h>
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

TString effLabel("PV");

double nBX= 1843714224.0;
double L = 3423000; //b^-1

void noiseSubtraction(const std::vector<bool> recreateDistributions = {false,false,false,false,false,false}, // {dataZB,dataEB,MC,EMs,EMd}
                      const std::vector<const char*> inputFile = {"dataZBInputFiles.txt","dataEBInputFiles.txt","dataMBInputFiles.txt","eposInputFiles.txt","emInputFiles_single.txt","emInputFiles_double.txt"},
                      const std::vector<const char*> sampleLabel = {"DATA_ZB","DATA_EB","DATA_MB","EPOSLHC","STARlight_single","STARlight_double"},
                      const std::vector<double> xSection = {-1.,-1.,-1.,5.61,10.0,0.041}, // cross section: {sigma_dummy,sigma_dummy,sigma_had,sigma_EM_single,sigma_EM_double} b
                      const char* colLabel = "XeXe 5.44 TeV",
                      const char* triggerName = "HLT_HIZeroBias_v1"
                 )
{
  const int nSamples = 6;
  if (recreateDistributions.size()!=nSamples) {cout << "[ERROR] recreateDistributions must have " << nSamples << " entries. Usage: {dataZB,dataEB,dataMB,MC,EMs,EMd}" << endl; return;}
  if (inputFile.size()!=nSamples) {cout << "[ERROR] inputfile must have " << nSamples << " entries. Usage: {dataZB,dataEB,dataMB,MC,EMs,EMd}" << endl; return;}
  if (sampleLabel.size()!=nSamples) {cout << "[ERROR] sampleLabel must have " << nSamples << " entries. Usage: {dataZB,dataEB,dataMB,MC,EMs,EMd}" << endl; return;}
  if (xSection.size()!=nSamples) {cout << "[ERROR] xSection must have " << nSamples << " entries. Usage: {dummy,dummy,dummy,MC,EMs,EMd}" << endl; return;}
  
  // Create directory to save the plots
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/Results";
  string histosDIR = mainDIR + "/../../DataHistos";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  TObjArray* arrHistosarr = new TObjArray();
  for (int n = 0 ; n < nSamples ; n++)
  {
    if (recreateDistributions[n])
    {
      cout << "[INFO] Creating distributions for " << sampleLabel[n]  << endl;
      mbDistributions(Form("%s/../../InputFiles/%s",gSystem->pwd(),inputFile[n]),triggerName,sampleLabel[n]);
    }
    
    TString tLabel("");
    if (!strcmp(sampleLabel[n],"DATA_ZB")) tLabel = "ZB";
    if (!strcmp(sampleLabel[n],"DATA_EB")) tLabel = "UB";
    if (!strcmp(sampleLabel[n],"DATA_MB") || !strcmp(sampleLabel[n],"EPOSLHC") || !strcmp(sampleLabel[n],"STARlight_single") || !strcmp(sampleLabel[n],"STARlight_double")) tLabel = "MB";
    
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
  
  const int nvar = 3;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFminus"};
  double varThr[nvar] = {35., 25., 25.};
  
  // Make the plots
  TObjArray* arrHistosZB = static_cast<TObjArray*>(arrHistosarr->At(0));
  TObjArray* arrHistosEB = static_cast<TObjArray*>(arrHistosarr->At(1));
  TObjArray* arrHistosMB = static_cast<TObjArray*>(arrHistosarr->At(2));
  TObjArray* arrHistosEMs = static_cast<TObjArray*>(arrHistosarr->At(4));
  TObjArray* arrHistosEMd = static_cast<TObjArray*>(arrHistosarr->At(5));
  
  TObjArray* arrVarZB(0x0);
  TObjArray* arrVarEB(0x0);
  TObjArray* arrVarMB(0x0);
  TObjArray* arrVarEMs(0x0);
  TObjArray* arrVarEMd(0x0);
  TH1* h_ZB(0x0);
  TH1* h_MB_sel(0x0);
  TH1* h_EB(0x0);
  TH1* h_EMs(0x0);
  TH1* h_EMd(0x0);
  for (int i = 0 ; i < nvar ; i++)
  {
    TString varN(varname[i]);
    
    
    if (arrHistosZB) arrVarZB = static_cast<TObjArray*>(arrHistosZB->FindObject(Form("histos_%s",varN.Data())));
    else {cout << "[ERROR] No ZB histos array for var " << varN.Data() << " was found" << endl; return;}
    if (arrHistosMB) arrVarMB = static_cast<TObjArray*>(arrHistosMB->FindObject(Form("histos_%s",varN.Data())));
    else {cout << "[ERROR] No MB histos array for var " << varN.Data() << " was found" << endl; return;}
    if (arrHistosEB) arrVarEB = static_cast<TObjArray*>(arrHistosEB->FindObject(Form("histos_%s",varN.Data())));
    else {cout << "[ERROR] No EB histos array for var " << varN.Data() << " was found" << endl; return;}
    if (arrHistosEMs) arrVarEMs = static_cast<TObjArray*>(arrHistosEMs->FindObject(Form("histos_%s",varN.Data())));
    else {cout << "[ERROR] No EMs histos array for var " << varN.Data() << " was found" << endl; return;}
    if (arrHistosEMd) arrVarEMd = static_cast<TObjArray*>(arrHistosEMd->FindObject(Form("histos_%s",varN.Data())));
    else {cout << "[ERROR] No EMd histos array for var " << varN.Data() << " was found" << endl; return;}
    
    h_ZB = static_cast<TH1*>(arrVarZB->FindObject(Form("h_ZB_%s",varN.Data()))->Clone());
    h_EB = static_cast<TH1*>(arrVarEB->FindObject(Form("h_UB_%s",varN.Data())));
    h_MB_sel = static_cast<TH1*>(arrVarMB->FindObject(Form("h_MB_BS_PV_HF3_%s",varN.Data()))->Clone());
    h_EMs = static_cast<TH1*>(arrVarEMs->FindObject(Form("h_%s",varN.Data()))->Clone());
    h_EMd = static_cast<TH1*>(arrVarEMd->FindObject(Form("h_%s",varN.Data()))->Clone());
    
    if (!h_ZB)
    {
      cout << "[ERROR] No ZB histos found in array " << endl;
      return;
    }
    if (!h_MB_sel)
    {
      cout << "[ERROR] No MB histo found in array " << endl;
      return;
    }
    if (!h_EMs || !h_EMd)
    {
      cout << "[ERROR] No EM histos found in array " << endl;
      return;
    }
    if (!h_EB) {cout << "[ERROR] No EB histo found in array " << endl; return;}
    
    // Normalisation factors
    double wZB = nBX/h_ZB->Integral();
    h_ZB->Scale(wZB);
    
    double wMB = h_ZB->Integral(h_ZB->FindBin(50.),h_ZB->GetNbinsX())/h_MB_sel->Integral(h_MB_sel->FindBin(50.),h_MB_sel->GetNbinsX());
    h_MB_sel->Scale(wMB);
    
    double sigma_MB = xSection[3] + xSection[4] + xSection[5];
    double nMBtrue = L*sigma_MB;
    double nNOISEtrue = nBX - nMBtrue;
    double wEB = nNOISEtrue/h_EB->Integral();
    h_EB->Scale(wEB);
    
    double wEMs = L*(xSection[4]/h_EMs->Integral());
    h_EMs->Scale(wEMs);
    
    double wEMd = L*(xSection[5]/h_EMd->Integral());
    h_EMd->Scale(wEMd);
    
    
    TH1* h_ZB_subs = static_cast<TH1*>(h_ZB->Clone("h_ZB_noiseSubs"));
    TH1* h_ZB_subs_emsubs = static_cast<TH1*>(h_ZB->Clone("h_ZB_noiseSubs_emSubs"));

    
    double xMin = 0.;
    double xMax = 50.;
    bool logY = true;
    if(!strcmp(varN.Data(),"trkEta") ) {xMin = -4.; xMax = 4.; logY = true;}
    if(!strcmp(varN.Data(),"trkPt") ) {xMin = 0.; xMax = 15.; logY = true;}
    
    //======= Draw Canvas 0 ((PV,NotPV),(combEff,effFactor))
        //** Pad 1
    gStyle -> SetOptStat(0);
    TCanvas* c0=  new TCanvas(Form("c0_%s",varN.Data()),"", 900,500);
    c0->Divide(2,1);
    c0->cd(1);
    if (logY) gPad->SetLogy();
    h_ZB->GetXaxis()->SetRangeUser(xMin,xMax);
    h_ZB->SetLineColor(1);
    h_ZB->SetMarkerColor(1);
    h_ZB->SetMarkerStyle(20);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      h_ZB->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      h_ZB->GetXaxis()->SetTitleOffset(1.3);
    }
    h_ZB->DrawClone();
    
    h_EB->SetLineColor(2);
    h_EB->SetMarkerColor(2);
    h_EB->SetMarkerStyle(21);
    h_EB->DrawClone("same");
    
    h_MB_sel->GetXaxis()->SetRangeUser(xMin,xMax);
    h_MB_sel->SetLineColor(3);
    h_MB_sel->SetMarkerColor(3);
    h_MB_sel->SetMarkerStyle(21);
    h_MB_sel->DrawClone("same");
    
    TLegend* l0 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l0);
    
    l0->AddEntry(h_ZB, "ZB", "p");
    l0->AddEntry(h_EB, "EB", "p");
    l0->AddEntry(h_MB_sel, "MB sel.", "p");
    l0->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(sampleLabel[0],0.54,0.2+0.60);
            //**
    
        //** Pad 2
    c0->cd(2);
    if (logY) gPad->SetLogy();
    h_ZB_subs->Add(h_EB,-1.);
    h_ZB_subs->GetXaxis()->SetRangeUser(xMin,xMax);
    h_ZB_subs->SetLineColor(1);
    h_ZB_subs->SetMarkerColor(1);
    h_ZB_subs->SetMarkerStyle(20);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      h_ZB_subs->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      h_ZB_subs->GetXaxis()->SetTitleOffset(1.3);
    }
    h_ZB_subs->Draw();
    //
    
    h_MB_sel->DrawClone("same");
    
    //
    
    h_EMs->SetLineColor(7);
    h_EMs->SetMarkerColor(7);
    h_EMs->SetMarkerStyle(26);
    h_EMs->Draw("same");
    h_EMd->SetLineColor(9);
    h_EMd->SetMarkerColor(9);
    h_EMd->SetMarkerStyle(32);
    h_EMd->Draw("same");

    
    TLegend* l01 = new TLegend(0.30,0.64,0.9,0.76);
    legStyle(l01);
    
    l01->AddEntry(h_ZB_subs, "ZB - EB", "p");
    l01->AddEntry(h_MB_sel, "MB sel.", "p");
    l01->AddEntry(h_EMs, "EM single", "p");
    l01->AddEntry(h_EMd, "EM double", "p");
    l01->Draw("same");
    
    c0->SaveAs(Form("%s/zbNoiseSubsNorm_%s.pdf",saveDIR.c_str(),varN.Data()));
    aSave->Add(c0);
        //**
    
//    c1->cd(1);
//
//    // Subtract EM events from ZB-Noise
//    h_ZB_subs_emsubs->Add(h_EMs,-1.);
//    h_ZB_subs_emsubs->Add(h_EMd,-1.);
//    
//    c1->cd(2);
//    h_ZB_subs_emsubs->Draw();
//    h_MB_sel->Draw("same");
//
//    
//    TLegend* l21 = new TLegend(0.40,0.64,1.00,0.76);
//    legStyle(l21);
//    l21->AddEntry(h_ZB_subs_emsubs,"ZB - Noise - EM", "p");
//    l21->AddEntry(h_MB_sel,"MB + BS + PV + HF3", "p");
//    
//    l21->Draw("same");
//    
//    drawText(Form("#varepsilon (selection + contam.) = %0.2f ",intMBTot/h_ZB_subs_emsubs->Integral(1,h_ZB_subs_emsubs->GetNbinsX())),0.45,0.20);
//    c1->SaveAs(Form("%s/zbNoiseSubs_%s.pdf",saveDIR.c_str(),varN.Data()));
//    aSave->Add(c1);
  }

  // save results
  TFile *fSave = new TFile(Form("%s/zbNoiseSubs_Norm.root",saveDIR.c_str()),"RECREATE");
  aSave->Write("zbNoiseSub", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
}
