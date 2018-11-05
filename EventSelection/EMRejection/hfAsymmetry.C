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

void hfAsymmetry(const std::vector<bool> recreateDistributions = {false,false,false}, // {had,EM}
                 const std::vector<const char*> inputFile = {"eposInputFiles.txt","emInputFiles_single.txt","dataMBInputFiles.txt"},
                 const std::vector<const char*> sampleLabel = {"EPOSLHC","STARlight_single","DATA_MB"},
                 const char* colLabel = "XeXe 5.44 TeV",
                 const char* triggerName = "HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1;HLT_HIL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1"
                 )
{
  
  TString tLabel(triggerName);
  if (tLabel.Contains("MinimumBias")) tLabel = "MB";
  else {cout << "[ERROR] Unkown trigger label" << endl; return;}
  
  // Check the input parameters
  const int nSamples = 3;
  if (recreateDistributions.size()!=nSamples) {cout << "[ERROR] recreateDistributions must have " << nSamples << " entries. Usage: {MC-had,starlight_single}" << endl; return;}
  if (inputFile.size()!=nSamples) {cout << "[ERROR] inputfile must have " << nSamples << " entries. Usage: {MC-had,starlight_single}" << endl; return;}
  if (sampleLabel.size()!=nSamples) {cout << "[ERROR] sampleLabel must have " << nSamples << " entries. Usage: {MC-had,starlight_single}" << endl; return;}
  
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
  
  // Create directory to save the plots
  string saveDIR = mainDIR + "/Results";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  const int nvar = 3;
  const char* varname[nvar] = {"hiHFasym","leadingTowerEAsym","leadingTowerETAsym"};
  
  const int nSel = 3;
  const char* selname[nSel] = {"", "MB", "MB_BS_PV_HF3"};
  
  // Define cross section weigthing functions
  TF1* fW = new TF1("wEMS","[0]",0.,1000.);
  
  // Make the plots
  TObjArray* arrVarHAD(0x0);
  TObjArray* arrVarEM(0x0);
  TObjArray* arrVarDATA(0x0);
  
  TH1* h_had(0x0);
  TH1* h_em(0x0);
  TH1* h_data(0x0);
  TH1* h_em_mirr(0x0);
  TH1* h_had_eff(0x0);
  TH1* h_em_eff(0x0);
    for (int k = 0 ; k < nvar ; k++)
    {
      TObjArray* arrHistosSave = new TObjArray();
      arrHistosSave->SetOwner(true);
      arrHistosSave->SetName(Form("histos_%s",varname[k]));
      
      TString xLabel("");
      TString xLabelAbs("");
      if (!strcmp("hiHFasym",varname[k]))
      {
        xLabel = "(#sum E_{T}^{+} - #sum E_{T}^{-}) / #sum E_{T}";
        xLabelAbs = "|#sum E_{T}^{+} - #sum E_{T}^{-}| / #sum E_{T}";
      }
      if (!strcmp("leadingTowerEAsym",varname[k]))
      {
        xLabel = "( E^{lT+} - E^{lT-}) / ( E^{lT+} + E^{lT-})";
        xLabelAbs = "|E^{lT+} - E^{lT-}| / ( E^{lT+} + E^{lT-})";
      }
      if (!strcmp("leadingTowerETAsym",varname[k]))
      {
        xLabel = "( E^{lT+}_{T} - E_{T}^{lT-}) / ( E^{lT+}_{T} + E_{T}^{lT-})";
        xLabelAbs = "|E^{lT+}_{T} - E_{T}^{lT-}| / ( E^{lT+}_{T} + E_{T}^{lT-})";
      }
      
      for (int i = 0 ; i < nSel ; i++)
      {
        TString selN(selname[i]);
        
        TObjArray* arrHistosHAD = static_cast<TObjArray*>(arrHistosarr->At(0));
        if (arrHistosHAD) arrVarHAD = static_cast<TObjArray*>(arrHistosHAD->FindObject(Form("histos_%s",varname[k])));
        else {cout << "[ERROR] No array of histos found for " << sampleLabel[0] << endl; return;}
        
        h_had = static_cast<TH1*>(arrVarHAD->FindObject(Form("h%s_%s",selN.IsNull()?"":Form("_%s",selN.Data()),varname[k]))->Clone());
        if (!h_had)
        {
          cout << "[ERROR] No hadronic histos found in array for selection " << selN.Data() << endl;
          return;
        }
        h_had->Sumw2();
        h_had_eff = static_cast<TH1*>(h_had->Clone(Form("%s_eff",h_had->GetName())));
        
        TObjArray* arrHistosEM = static_cast<TObjArray*>(arrHistosarr->At(1));
        if (arrHistosEM) arrVarEM = static_cast<TObjArray*>(arrHistosEM->FindObject(Form("histos_%s",varname[k])));
        else {cout << "[ERROR] No array of histos found for " << sampleLabel[1] << endl; return;}
        
        h_em = static_cast<TH1*>(arrVarEM->FindObject(Form("h%s_%s",selN.IsNull()?"":Form("_%s",selN.Data()),varname[k]))->Clone());
        if (!h_em)
        {
          cout << "[ERROR] No em histos found in array for selection " << selN.Data() << endl;
          return;
        }
        h_em->Sumw2();
        h_em_eff = static_cast<TH1*>(h_em->Clone(Form("%s_eff",h_em->GetName())));
        h_em_mirr = static_cast<TH1*>(h_em->Clone(Form("%s_mirr",h_em->GetName())));
        
        TObjArray* arrHistosData = static_cast<TObjArray*>(arrHistosarr->At(2));
        if (arrHistosData) arrVarDATA = static_cast<TObjArray*>(arrHistosData->FindObject(Form("histos_%s",varname[k])));
        else {cout << "[ERROR] No array of histos found for " << sampleLabel[2] << endl; return;}
        
        h_data = static_cast<TH1*>(arrVarDATA->FindObject(Form("h%s_%s",selN.IsNull()?"":Form("_%s",selN.Data()),varname[k]))->Clone());
        if (!h_data)
        {
          cout << "[ERROR] No em histos found in array for selection " << selN.Data() << endl;
          return;
        }
        h_data->Sumw2();
        
        int nBins = h_had->GetNbinsX();
        
        // Make full EM distribution by adding the mirror distribution
        for (int j = 1 ; j<nBins ; j++)
        {
          double binC = h_em->GetBinCenter(j);
          int binPlus = h_em->FindBin(-binC);
          
          h_em_mirr->SetBinContent(binPlus,h_em->GetBinContent(j));
          h_em_mirr->SetBinError(binPlus,h_em->GetBinError(j));
        }
        h_em_mirr->Add(h_em);
        
        gStyle -> SetOptStat(0);
        TCanvas* c=  new TCanvas(Form("c%s_%s",selN.IsNull()?"":Form("_%s",selN.Data()),varname[k]),"", 900,500);
        c->Divide(2,1);
        
        c->cd(1);
        gPad->SetLogy();
        h_had->GetXaxis()->SetRangeUser(-1.,1.);
        h_had->Scale(1./h_had->Integral(1,nBins));
        h_had->SetLineColor(1);
        h_had->SetMarkerColor(1);
        h_had->SetMarkerStyle(20);
        h_had->SetXTitle(xLabel.Data());
        h_had->GetXaxis()->SetTitleOffset(1.3);
        h_had->DrawClone();
        
        h_data->Scale(1./h_data->Integral(1,nBins));
        h_data->SetLineColor(4);
        h_data->SetMarkerColor(4);
        h_data->SetMarkerStyle(21);
        h_data->DrawClone("same");
        
        h_em->Scale(1./h_em->Integral(1,nBins));
        h_em->SetLineColor(2);
        h_em->SetMarkerColor(2);
        h_em->SetMarkerStyle(20);
        h_em->DrawClone("same");
        
        h_em_mirr->Scale(1./h_em_mirr->Integral(1,nBins));
        h_em_mirr->SetLineColor(2);
        h_em_mirr->SetMarkerColor(2);
        h_em_mirr->SetMarkerStyle(24);
        h_em_mirr->DrawClone("same");
        
        TLegend* l1 = new TLegend(0.25,0.25,0.85,0.37);
        legStyle(l1);
        l1->AddEntry(h_had, sampleLabel[0], "p");
        l1->AddEntry(h_em, sampleLabel[1], "p");
        l1->AddEntry(h_em_mirr, Form("%s_mirr",sampleLabel[1]), "p");
        l1->AddEntry(h_data, sampleLabel[2], "p");
        l1->Draw("same");
        
        drawText(colLabel,0.54,0.2+0.65);
        drawText(Form("%s filters",selN.IsNull()?"No":selN.Data()),0.54,0.2+0.60);
        
        // Compute efficiency vs hfAsym
        double intHad = h_had->Integral(1,nBins);
        double intEM = h_em_mirr->Integral(1,nBins);
        
        for (int j = 1 ; j<nBins ; j++)
        {
          double binC = h_had->GetBinCenter(j);
          int binPlus = h_had->FindBin(-binC);
          
          double eff_had = h_had->Integral(j,binPlus)/intHad; // Efficiency in keeping had events
          double eff_em = (h_em_mirr->Integral(1,j) + h_em_mirr->Integral(binPlus,nBins))/intEM;  // Efficiency in rejecting EM events
          
          if (binC>0)
          {
            h_had_eff->SetBinContent(binPlus,0);
            h_had_eff->SetBinError(binPlus,0);
            //        cout << binC << " ; " << binPlus << " ; " << eff_em << endl;
            h_em_eff->SetBinContent(binPlus,0);
            h_em_eff->SetBinError(binPlus,0);
          }
          else{
            h_had_eff->SetBinContent(binPlus,eff_had);
            h_had_eff->SetBinError(binPlus,0);
            //        cout << binC << " ; " << binPlus << " ; " << eff_em << endl;
            h_em_eff->SetBinContent(binPlus,eff_em);
            h_em_eff->SetBinError(binPlus,0);
          }
        }
        
        c->cd(2);
        h_had_eff->GetXaxis()->SetRangeUser(0.,1.);
        h_had_eff->GetYaxis()->SetRangeUser(0.,1.1);
        h_had_eff->SetLineColor(1);
        h_had_eff->SetMarkerColor(1);
        h_had_eff->SetMarkerStyle(20);
        h_had_eff->SetXTitle(xLabelAbs.Data());
        h_had_eff->GetXaxis()->SetTitleOffset(1.3);
        h_had_eff->DrawClone();
        
        h_em_eff->GetXaxis()->SetRangeUser(0.,1.);
        h_em_eff->SetLineColor(2);
        h_em_eff->SetMarkerColor(2);
        h_em_eff->SetMarkerStyle(20);
        h_em_eff->DrawClone("same");
        
        TLegend* l2 = new TLegend(0.35,0.64,0.95,0.76);
        legStyle(l2);
        l2->AddEntry(h_had_eff, Form("%s #varepsilon_{cut}",sampleLabel[0]), "l");
        l2->AddEntry(h_em_eff, Form("%s 1 - #varepsilon_{cut}",sampleLabel[1]), "l");
        l2->Draw("same");
        
        c->SaveAs(Form("%s/hfAsymRejection_%s_%s.pdf",saveDIR.c_str(),selN.IsNull()?"NoSel":selN.Data(),varname[k]));
        arrHistosSave->Add(c);
      }
      aSave->Add(arrHistosSave);
    }

  // save results
  TFile *fSave = new TFile(Form("%s/hfAsymRejection.root",saveDIR.c_str()),"RECREATE");
  aSave->Write("hfAsymRejection", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close();
}
