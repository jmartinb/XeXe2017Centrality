#ifndef treeUtils
#define treeUtils
#include <TFile.h>
#include <TTree.h>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCollection.h"
#include "TDirectoryFile.h"

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;

//////////////////////////////////////////////
// Definition of global variables
//////////////////////////////////////////////

float hf, hfplus, hfpluseta4, hfminuseta4, hfminus, hfhit,zvtx[10000], trkdxy1[10000], trkdz1[10000],normchi2vtx[10000], vtxdist2d[10000],trkEta[10000],trkPt[10000];
int nTowers,hiBin,hiNtrks,ntrk,npix,nvtx,maxmultvtx,maxptvtx,ntrkvtx[10000];
float towerE[4000], towerET[4000], eta[4000];
int ieta[4000];
int ProcessID;
bool highpurity[10000];
UInt_t run, lumi;
int HLT_Trigger_part1, HLT_Trigger_part2, HLT_Trigger_part3, HLT_Trigger_part4, HLT_Trigger_part5, HLT_Trigger_part6, HLT_Trigger_part7, HLT_Trigger_part8, HLT_Trigger_part9, HLT_Trigger_part10, HLT_Trigger_part11, HLT_Trigger_part12, HLT_Trigger_part13, HLT_Trigger_part14, HLT_Trigger_part15, HLT_Trigger_part16, HLT_Trigger_part17, HLT_Trigger_part18, HLT_Trigger_part19, HLT_Trigger_part20;
int pBeamScrapingFilter, pclusterCompatibilityFilter, pprimaryVertexFilter, phfCoincFilter1, phfCoincFilter2, phfCoincFilter3, phfCoincFilter4, phfCoincFilter5;


//////////////////////////////////////////////
//////////////////////////////////////////////

TChain* readFiles(const char* inputFile)
{
  // Read files and chain TTrees
  TChain *tref;
  char line[1024];
  ifstream in(inputFile);
  tref = new TChain("hiEvtAnalyzer/HiTree","");
  TChain * tskimanalysis = new TChain("skimanalysis/HltTree","");
  TChain * thltanalysis = new TChain("hltanalysisReco/HltTree",""); // This name is specific for XeXe, use hltanalysis otherwise
  TChain * ttrack= new TChain("ppTrack/trackTree","");
  TChain * trechits= new TChain("rechitanalyzer/tower","");
  int cnt(0);
  bool eTrackTree(kFALSE);
  bool eRechitsTree(kFALSE);
  while (in.getline(line,1024,'\n'))
  {
    if (cnt==0) // Check if tracks and rechits are there
    {
      TFile* t = TFile::Open(line,"READ");
      TDirectoryFile* dirTrk = static_cast<TDirectoryFile*>(t->FindObjectAny("ppTrack"));
      if (dirTrk) eTrackTree = kTRUE;
      TDirectoryFile* dirRH = static_cast<TDirectoryFile*>(t->FindObjectAny("rechitanalyzer"));
      if (dirRH) eRechitsTree = kTRUE;
      t->Close();
    }
    
    const char* fileName = Form("%s",line);
    tref->Add(fileName);
    tskimanalysis->Add(fileName);
    thltanalysis->Add(fileName);
    if (eTrackTree) ttrack->Add(fileName);
    if (eRechitsTree) trechits->Add(fileName);
    cnt++;
  }
  
  tref->AddFriend(tskimanalysis);
  tref->AddFriend(thltanalysis);
  if (eTrackTree) tref->AddFriend(ttrack);
  if (eRechitsTree) tref->AddFriend(trechits);
  
  return tref;
}
//////////////////////////////////////////////
//////////////////////////////////////////////

bool loadBranches(const char* triggerName, TChain* t)
{
  TString tName(triggerName);
  
  TObjArray* tArray(0x0);
  if (tName.Contains(";")) tArray = tName.Tokenize(";");
  else
  {
    tArray = new TObjArray();
    tArray->Add(new TObjString(triggerName));
  }
  tArray->SetOwner(kTRUE);
  
  TIter it(tArray);
  TObjString* obj(0x0);
  int cnt(1);
  while ((obj = static_cast<TObjString*>(it.Next())))
  {
    TString tChar = obj->GetString();
    if (!t->GetBranch(tChar.Data()))
    {
      cout << "[ERROR] Specified trigger path: " << tChar.Data() << " does not exist" << endl;
      return false;
    }
    else if (cnt==1) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part1);
    else if (cnt==2) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part2);
    else if (cnt==3) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part3);
    else if (cnt==4) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part4);
    else if (cnt==5) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part5);
    else if (cnt==6) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part6);
    else if (cnt==7) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part7);
    else if (cnt==8) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part8);
    else if (cnt==9) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part9);
    else if (cnt==10) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part10);
    else if (cnt==11) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part11);
    else if (cnt==12) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part12);
    else if (cnt==13) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part13);
    else if (cnt==14) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part14);
    else if (cnt==15) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part15);
    else if (cnt==16) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part16);
    else if (cnt==17) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part17);
    else if (cnt==18) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part18);
    else if (cnt==19) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part19);
    else if (cnt==20) t->SetBranchAddress(tChar.Data(),	&HLT_Trigger_part20);
    else { cout << "[ERROR] Code does not support such a long trigger list, please update" << endl; return false;}
    
    cnt++;
  }
  
  t->SetBranchAddress("run",	&run);
  t->SetBranchAddress("lumi",	&lumi);
  if (t->GetBranch("ProcessID")) t->SetBranchAddress("ProcessID",	&ProcessID);
  t->SetBranchAddress("hiHF",		&hf);
  t->SetBranchAddress("hiBin",		&hiBin);
  t->SetBranchAddress("hiHFplus",	&hfplus);
  t->SetBranchAddress("hiHFplusEta4",	&hfpluseta4);
  t->SetBranchAddress("hiHFminus",	&hfminus);
  t->SetBranchAddress("hiHFminusEta4",	&hfminuseta4);
  t->SetBranchAddress("hiNtracks",	&hiNtrks);
  t->SetBranchAddress("hiHFhit",	&hfhit);
  t->SetBranchAddress("hiNpix",		&npix);
  if (t->GetBranch("nVtx")) t->SetBranchAddress("nVtx",		&nvtx);
  if (t->GetBranch("nTrk")) t->SetBranchAddress("nTrk",		&ntrk);
  if (t->GetBranch("trkEta")) t->SetBranchAddress("trkEta",		&trkEta);
  if (t->GetBranch("trkPt")) t->SetBranchAddress("trkPt",		&trkPt);
  if (t->GetBranch("maxMultVtx")) t->SetBranchAddress("maxMultVtx",		&maxmultvtx);
  if (t->GetBranch("maxPtVtx")) t->SetBranchAddress("maxPtVtx",		&maxptvtx);
  if (t->GetBranch("zVtx")) t->SetBranchAddress("zVtx",		&zvtx);
  if (t->GetBranch("trkDxy1")) t->SetBranchAddress("trkDxy1",		&trkdxy1);
  if (t->GetBranch("trkDz1")) t->SetBranchAddress("trkDz1",		&trkdz1);
  if (t->GetBranch("nTrkVtx")) t->SetBranchAddress("nTrkVtx",		&ntrkvtx);
  if (t->GetBranch("highPurity")) t->SetBranchAddress("highPurity",		&highpurity);
  if (t->GetBranch("normChi2Vtx")) t->SetBranchAddress("normChi2Vtx",		&normchi2vtx);
  if (t->GetBranch("vtxDist2D")) t->SetBranchAddress("vtxDist2D",		&vtxdist2d);
  if (t->GetBranch("pBeamScrapingFilter")) t->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);
  if (t->GetBranch("pclusterCompatibilityFilter"))  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  if (t->GetBranch("pPAprimaryVertexFilter")) t->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter);
  if (t->GetBranch("pprimaryVertexFilter")) t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  if (t->GetBranch("phfCoincFilter1")) t->SetBranchAddress("phfCoincFilter1", &phfCoincFilter1);
  if (t->GetBranch("phfCoincFilter2")) t->SetBranchAddress("phfCoincFilter2", &phfCoincFilter2);
  if (t->GetBranch("phfCoincFilter3")) t->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3);
  if (t->GetBranch("phfCoincFilter4")) t->SetBranchAddress("phfCoincFilter4", &phfCoincFilter4);
  if (t->GetBranch("phfCoincFilter5")) t->SetBranchAddress("phfCoincFilter5", &phfCoincFilter5);
  
  if (t->GetBranch("n")) t->SetBranchAddress("n", &nTowers);
  if (t->GetBranch("e")) t->SetBranchAddress("e", &towerE);
  if (t->GetBranch("et")) t->SetBranchAddress("et", &towerET);
  if (t->GetBranch("ieta")) t->SetBranchAddress("ieta", &ieta);
  if (t->GetBranch("eta")) t->SetBranchAddress("eta", &eta);
  
  tArray->Delete();
  return true;
}

void MergeTrees(TChain& tree, std::ifstream& filelist)
{
  if (not (filelist.is_open())) {
    std::cerr << "failed in open." << std::endl;
    exit(EXIT_FAILURE);
  }

  char line[BUFSIZ];
  int counter = 0; 

  while (not (filelist.eof())) {
    filelist.getline(line, sizeof(line));

    if (strcmp(line, "") == 0) continue;

    TFile a(line);
    if (!a.IsOpen()) {
      std::cerr << "failed in open :" << line << std::endl;
      exit(EXIT_FAILURE);
    }
    if (gROOT->FindObject(tree.GetName()) == 0) {
      std::cerr << tree.GetName() << " is not found in :" << line << std::endl;
      exit(EXIT_FAILURE);
    }

    counter ++ ;
    tree.AddFile(line);
  }

  std::cout << counter << " files are loaded." << std::endl;
  return;
}

#endif
