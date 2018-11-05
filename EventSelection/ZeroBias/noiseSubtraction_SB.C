// Subtraction of noise from ZB sample
// This macro uses a sideband method to create an histogram of the noise from ZB
// Author : Javier Martin Blanco 02/03/2018

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

void noiseSubtraction_SB(const std::vector<bool> recreateDistributions = {false,false,false,false,false,false}, // {dataZB,dataEB,MC,EMs,EMd}
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
    if (!strcmp(sampleLabel[n],"DATA_EB")) tLabel = "NBOR";
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
  
  const int nvar = 5;
  const char* varname[nvar] = {"hiHF", "hiHFplus", "hiHFplusEta4", "hiHFminus","hiHFminusEta4"};
  double varThr[nvar] = {35., 25., 20., 25., 20.};
  
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
  TH1* h_PV(0x0);
  TH1* h_NotPV(0x0);
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
    h_PV = static_cast<TH1*>(arrVarZB->FindObject(Form("h_ZB_BS_PV_%s",varN.Data()))->Clone());
    h_NotPV = static_cast<TH1*>(arrVarZB->FindObject(Form("h_ZB_BS_NotPV_%s",varN.Data()))->Clone());
    h_MB_sel = static_cast<TH1*>(arrVarMB->FindObject(Form("h_MB_BS_PV_HF3_%s",varN.Data()))->Clone());
    h_EMs = static_cast<TH1*>(arrVarEMs->FindObject(Form("h_%s",varN.Data()))->Clone());
    h_EMd = static_cast<TH1*>(arrVarEMd->FindObject(Form("h_%s",varN.Data()))->Clone());
    
    if (!h_ZB || !h_PV || !h_NotPV)
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
    
    
    TH1* h_noise = static_cast<TH1*>(h_NotPV->Clone("h_noise"));
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
    
    h_PV->SetLineColor(2);
    h_PV->SetMarkerColor(2);
    h_PV->SetMarkerStyle(21);
    h_PV->DrawClone("same");
    
    h_NotPV->GetXaxis()->SetRangeUser(xMin,xMax);
    h_NotPV->SetLineColor(3);
    h_NotPV->SetMarkerColor(3);
    h_NotPV->SetMarkerStyle(21);
    h_NotPV->DrawClone("same");
    
    TLegend* l0 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l0);
    
    l0->AddEntry(h_ZB, "ZB", "p");
    l0->AddEntry(h_PV, "ZB + PV", "p");
    l0->AddEntry(h_NotPV, "ZB + not PV", "p");
    l0->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(sampleLabel[0],0.54,0.2+0.60);
            //**
    
    
    // Create combined PV eff (had + EM) and histo (1-eff)/eff
    TEfficiency* eff_EPOS = getMCeff(varN.Data(),sampleLabel[3],"0","",effLabel.Data(),"");
    TEfficiency* eff_StarlightS = getMCeff(varN.Data(),sampleLabel[4],"0","",effLabel.Data(),"");
    TEfficiency* eff_StarlightD = getMCeff(varN.Data(),sampleLabel[5],"0","",effLabel.Data(),"");
    
    if (!eff_EPOS || !eff_StarlightS || !eff_StarlightD) {cout << "[ERROR] No efficiencies found for var " << varN.Data() << endl; return;}
    
    double sigma_MB = xSection[3] + xSection[4] + xSection[5];
    TH1* eff_factor = eff_EPOS->GetCopyPassedHisto();
    TH1* eff_comb = eff_EPOS->GetCopyPassedHisto();
    for (int bin = 1 ; bin < eff_comb->GetNbinsX() ; bin++)
    {
      double eff_had = eff_EPOS->GetEfficiency(bin);
      double eff_emS = eff_StarlightS->GetEfficiency(bin);
      double eff_emD = eff_StarlightD->GetEfficiency(bin);
      
      double eff = (xSection[3]*eff_had + xSection[4]*eff_emS + xSection[5]*eff_emD)/sigma_MB;
      if (eff<=0.) // In this case we are at very low energy and all the efficiencies are 0
      {
        eff_factor->SetBinContent(bin,0.);
        eff_factor->SetBinError(bin,0.);
        eff_comb->SetBinContent(bin,0.);
        eff_comb->SetBinError(bin,0.);
        continue;
      }
      
      if (eff_emS == 0 && eff_emD == 0) // To avoid problems when the EM efficiencies cannot be calculated at high energy (no Starlight statistics)
      {
        eff_emS = 1.0;
        eff_emD = 1.0;
      }
      if (eff_emS == 0 && eff_emD != 0) eff_emS = 1.0;
      if (eff_emS != 0 && eff_emD == 0) eff_emD = 1.0;
      
      eff = (xSection[3]*eff_had + xSection[4]*eff_emS + xSection[5]*eff_emD)/sigma_MB;
      
      eff_factor->SetBinContent(bin,(1.-eff)/eff);
      eff_factor->SetBinError(bin,0.);
      eff_comb->SetBinContent(bin,eff);
      eff_comb->SetBinError(bin,0);
    }
    
        //** Pad 2
    c0->cd(2);
    eff_factor->GetXaxis()->SetRangeUser(xMin,xMax);
    eff_factor->GetYaxis()->SetRangeUser(0,2.5);
    eff_factor->SetLineColor(1);
    eff_factor->SetMarkerColor(1);
    eff_factor->SetMarkerStyle(25);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      eff_factor->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      eff_factor->GetXaxis()->SetTitleOffset(1.3);
    }
    eff_factor->DrawClone("P");
    
    eff_comb->SetLineColor(1);
    eff_comb->SetMarkerColor(1);
    eff_comb->SetMarkerStyle(24);
    eff_comb->DrawClone("sameP");
    
    TLegend* l01 = new TLegend(0.30,0.64,0.9,0.76);
    legStyle(l01);
    
    l01->AddEntry(eff_comb, "#varepsilon_{had+EM}", "p");
    l01->AddEntry(eff_factor, "(1-#varepsilon_{had+EM})/#varepsilon_{had+EM}", "p");
    l01->Draw("same");
    
    drawText("EPOS-LHC & STARlight",0.54,0.2+0.65);
    
    c0->SaveAs(Form("%s/cutDistrANDeff_%s.pdf",saveDIR.c_str(),varN.Data()));
    aSave->Add(c0);
        //**
    
    // Obtain noise distribution
    h_PV->Multiply(eff_factor);
    h_noise->Add(h_PV,-1.);
    
    //======= Draw Canvas ((ZB,Noise),((Noise,EB),ratio))
        //** Pad 1
    TCanvas* c=  new TCanvas(Form("c_%s",varN.Data()),"", 900,500);
    c->Divide(2,1);
    c->cd(1);
    gPad->SetTopMargin(0.09);
    if (logY) gPad->SetLogy();
    h_ZB->GetXaxis()->SetRangeUser(xMin,xMax);
    h_ZB->SetLineColor(1);
    h_ZB->SetMarkerColor(1);
    h_ZB->SetMarkerStyle(21);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      h_ZB->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      h_ZB->GetXaxis()->SetTitleOffset(1.3);
    }
    h_ZB->DrawCopy();

    h_noise->SetLineColor(2);
    h_noise->SetMarkerColor(2);
    h_noise->SetMarkerStyle(20);
    h_noise->DrawCopy("same");
    
    TLegend* l1 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l1);
    
    l1->AddEntry(h_ZB, "ZB", "p");
    l1->AddEntry(h_noise, "Noise (ZB_{NotPV} - #frac{1-#varepsilon_{PV}}{#varepsilon_{PV}}*ZB_{PV})", "p");
    l1->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(sampleLabel[0],0.54,0.2+0.60);
        //**
    
        //** Pad 2
    c->cd(2)->Divide(1,2);
    c->cd(2)->cd(1);
    gPad->SetPad(0., 0.3, 0.99, 0.99);
    gPad->SetTopMargin(0.09);
    if (logY) gPad->SetLogy();
    TH1* h_noise_norm = static_cast<TH1*>(h_noise->Clone("noise_norm"));
    h_noise_norm->Scale(1./h_noise_norm->Integral(1,h_noise_norm->GetNbinsX()));
    h_noise_norm->GetXaxis()->SetRangeUser(xMin,xMax);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      h_noise_norm->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      h_noise_norm->GetXaxis()->SetTitleOffset(1.3);
    }
    h_noise_norm->DrawCopy();
    h_EB = static_cast<TH1*>(arrVarEB->FindObject(Form("h_NBOR_%s",varN.Data())));
    if (!h_EB) {cout << "[ERROR] No EB histo found in array " << endl; return;}
    h_EB->SetLineColor(6);
    h_EB->SetMarkerColor(6);
    h_EB->SetMarkerStyle(20);
    h_EB->Scale(1./h_EB->Integral(1,h_EB->GetNbinsX()));
    h_EB->DrawCopy("same");
    
    TLegend* l12 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l12);
    l12->AddEntry(h_noise_norm, "Noise", "p");
    l12->AddEntry(h_EB, "EB", "p");
    l12->Draw("same");
    
    c->cd(2)->cd(2);
    gPad->SetPad(0., 0., 0.99, 0.3);
    gPad->SetBottomMargin(0.1/0.3);
    gPad->SetGridy();
    TH1* h_noise_ratio = static_cast<TH1*>(h_noise_norm->Clone("noise_ratio"));
    h_noise_ratio->Divide(h_EB);
    h_noise_ratio->GetXaxis()->SetRangeUser(xMin,xMax);
    h_noise_ratio->GetYaxis()->SetRangeUser(0.2,1.8);
    h_noise_ratio->GetXaxis()->SetLabelSize(0.09);
    h_noise_ratio->GetXaxis()->SetTitleSize(0.09);
    h_noise_ratio->GetYaxis()->SetLabelSize(0.09);
    h_noise_ratio->GetYaxis()->SetTitleSize(0.09);
    h_noise_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_noise_ratio->SetYTitle("Noise/EB");
    h_noise_ratio->SetLineColor(6);
    h_noise_ratio->SetMarkerColor(6);
    h_noise_ratio->SetMarkerStyle(20);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      h_noise_ratio->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      h_noise_ratio->GetXaxis()->SetTitleOffset(1.3);
    }
    h_noise_ratio->Draw();
    
    c->SaveAs(Form("%s/zbVSebNoise_%s.pdf",saveDIR.c_str(),varN.Data()));
    aSave->Add(c);
        //**
    
    
    //======= Draw Canvas 1 ((ZB-Noise,EMs,EMD),(ZB-Noise - EMs - EMD,MB selected))
    //** Pad 1
    TCanvas* c1=  new TCanvas(Form("c1_%s",varN.Data()),"", 900,500);
    c1->Divide(2,1);
    c1->cd(1);
    if (logY) gPad->SetLogy();
    h_ZB_subs->Add(h_noise,-1.);
    h_ZB_subs->GetXaxis()->SetRangeUser(xMin,xMax);
    h_ZB_subs->SetLineColor(3);
    h_ZB_subs->SetMarkerColor(3);
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
    TLegend* l2 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l2);
    
    l2->AddEntry(h_ZB_subs, "ZB - Noise", "p");
    l2->AddEntry(h_EMs, "EM single", "p");
    l2->AddEntry(h_EMd, "EM double", "p");
    l2->Draw("same");
    drawText(colLabel,0.54,0.2+0.65);
    drawText(sampleLabel[0],0.54,0.2+0.60);
    //
    
    c1->cd(2);
    if (logY) gPad->SetLogy();
    h_ZB_subs_emsubs->Add(h_noise,-1.);
    
    double intMB = h_MB_sel->Integral(h_MB_sel->FindBin(varThr[i]),h_MB_sel->GetNbinsX());
    double intZBsubs = h_ZB_subs_emsubs->Integral(h_ZB_subs_emsubs->FindBin(varThr[i]),h_ZB_subs_emsubs->GetNbinsX());
    
    h_ZB_subs_emsubs->GetXaxis()->SetRangeUser(xMin,xMax);
    h_ZB_subs_emsubs->SetLineColor(3);
    h_ZB_subs_emsubs->SetMarkerColor(3);
    h_ZB_subs_emsubs->SetMarkerStyle(20);
    if(varN.Contains("hiHF") && !varN.Contains("hit"))
    {
      TString axis(varN.Data());
      axis.Replace(0,2,"");
      h_ZB_subs_emsubs->SetXTitle(Form("%s #sum E_{T} (GeV)",axis.Data()));
      h_ZB_subs_emsubs->GetXaxis()->SetTitleOffset(1.3);
    }
   
    
    
    // Normalise MB to ZB
    h_MB_sel->Scale(intZBsubs/intMB);
    h_MB_sel->GetXaxis()->SetRangeUser(xMin,xMax);
    h_MB_sel->SetLineColor(4);
    h_MB_sel->SetMarkerColor(4);
    h_MB_sel->SetMarkerStyle(20);

    
    // Normalise EM to MB
    double intMBTot = h_ZB_subs_emsubs->Integral(1,h_ZB_subs_emsubs->GetNbinsX()); // NOTE: MB is actually ZB - Noise (so all had+EM events)
    double normEMs = (intMBTot/sigma_MB)*(xSection[4]/h_EMs->Integral(1,h_EMs->GetNbinsX()));
    h_EMs->Scale(normEMs);
    double normEMd = (intMBTot/sigma_MB)*(xSection[5]/h_EMd->Integral(1,h_EMd->GetNbinsX()));
    h_EMd->Scale(normEMd);
    
    c1->cd(1);
    h_EMs->SetLineColor(7);
    h_EMs->SetMarkerColor(7);
    h_EMs->SetMarkerStyle(26);
    h_EMs->Draw("same");
    h_EMd->SetLineColor(9);
    h_EMd->SetMarkerColor(9);
    h_EMd->SetMarkerStyle(32);
    h_EMd->Draw("same");
    
    // Subtract EM events from ZB-Noise
    h_ZB_subs_emsubs->Add(h_EMs,-1.);
    h_ZB_subs_emsubs->Add(h_EMd,-1.);
    
    c1->cd(2);
    h_ZB_subs_emsubs->Draw();
    h_MB_sel->Draw("same");

    
    TLegend* l21 = new TLegend(0.40,0.64,1.00,0.76);
    legStyle(l21);
    l21->AddEntry(h_ZB_subs_emsubs,"ZB - Noise - EM", "p");
    l21->AddEntry(h_MB_sel,"MB + BS + PV + HF3", "p");
    
    l21->Draw("same");
    
    drawText(Form("#varepsilon (selection + contam.) = %0.2f ",intMBTot/h_ZB_subs_emsubs->Integral(1,h_ZB_subs_emsubs->GetNbinsX())),0.45,0.20);
    c1->SaveAs(Form("%s/zbNoiseSubs_%s.pdf",saveDIR.c_str(),varN.Data()));
    aSave->Add(c1);
  }

  // save results
  TFile *fSave = new TFile(Form("%s/zbNoiseSubs.root",saveDIR.c_str()),"RECREATE");
  aSave->Write("zbNoiseSub", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
}
