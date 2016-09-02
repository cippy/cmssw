#define EoverP_cxx
#include "EoverP.h"

#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                           
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file                          
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators 

#include <Rtypes.h> // to use kColor

#include "histoFunc.h"
#include "Option.h"


#define DATA2016 1
#define SKIM_1LEP1JET_80X 1
#define FIT_2SIDE_CB 0      // 0 for single tail Crystal Ball for the fit, 1 for double tail
#define SET_SCALE_ON_Y 0    // select default or user defined ranges for y axis
#define USE_E 0 // 0 for ET and 1 for E in the binning

using namespace std;


//=====================================================================

// class Option {

// public:
//   Option();
//   ~Option();
//   Int_t Get_data2016() const { return data2016;}
//   Int_t Get_skim1lep1jet80X() const { return skim1lep1jet80X;}
//   Int_t Get_fit2sideCB() const { return fit2sideCB;}
//   Int_t Get_setScaleOnY() const { return setScaleOnY;}
//   Int_t Get_useE() const { return useE;}
//   string Get_dirName() const { return dirName;}

//   void Set_data2016(const Int_t &value) { data2016 = value; }
//   void Set_skim1lep1jet80X(const Int_t &value) { skim1lep1jet80X = value; }
//   void Set_fit2sideCB(const Int_t &value) { fit2sideCB = value; }
//   void Set_setScaleOnY(const Int_t &value) { setScaleOnY = value; }
//   void Set_useE(const Int_t &value) { useE = value; }
//   void Set_dirName(const string &value) { dirName = value; }


// private:
//   Int_t data2016;
//   Int_t skim1lep1jet80X;
//   Int_t fit2sideCB;
//   Int_t setScaleOnY;
//   Int_t useE;
//   string dirName;

// };  

// Option::Option() {
//   data2016 = 1;
//   skim1lep1jet80X = 1;
//   fit2sideCB = 0;
//   setScaleOnY = 0;
//   useE = 0;
//   dirName = "default";
// }

// Option::~Option() {
//   cout<<"Option::~Option() called"<<endl;
// }

//=====================================================================

void getEoP_templateMC(TH1F *hinput = NULL, const string& dirName = "") {

  string fileName = dirName + "EoverP_WJetsToLNu.root";  // template is created when running on MC, and is stored in this file

  TH1::AddDirectory(kFALSE);  // with this line the histogram taken with TFile::Get() is no more associated to the file (no need to Clone anymore)
  // note that even though I clone the object, if the file is closed the clone is lost, as it happens for the object retrieved with TFile::Get()

  TH1F *htmp = NULL;

  cout << "Opening file " << fileName << endl;

  TFile* f = TFile::Open(fileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  // get E/P template from file
  htmp = (TH1F*)f->Get("hEoP_template");  
  if (!htmp) {
    cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileName << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  htmp->SetDirectory(0);
  //  hinput = (TH1F*) htmp->Clone();
  hinput = new TH1F(*htmp);  // it looks like that after exiting the function hinput is gone
  if (!hinput) {
    cout << "Error in function getEoP_templateMC(TH1F * hinput, const string& dirName): hinput is NULL. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  //  f->Close();

}


//=====================================================================


void getLowAndUpBinGivenRange(const TH1F* histo, const Double_t &lowerX, const Double_t &upperX, Int_t &lowBin, Int_t &upBin) {
  lowBin = histo->GetXaxis()->FindFixBin(lowerX);
  upBin = histo->GetXaxis()->FindFixBin(upperX);
}

Int_t getLowBinGivenRange(const TH1F* histo, const Double_t &lowerX) {
  Int_t lowBin = histo->GetXaxis()->FindFixBin(lowerX);
  return lowBin;
}

Int_t getUpBinGivenRange(const TH1F* histo, const Double_t &upperX) {
  Int_t upBin = histo->GetXaxis()->FindFixBin(upperX);
  return upBin;
}

//=====================================================================

Double_t myCrystalBallRightTail(double* x, double* par) {

  Double_t xcur = x[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t invAbsAlpha = 1./absAlpha;

  if ( t <= absAlpha)  {   // would be t >= -absAlpha for left tail
    return N*exp(-0.5*t*t);
  } else {
    Double_t A = TMath::Power(n*invAbsAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t B = n*invAbsAlpha - absAlpha;
    return N*A * TMath::Power(B+t,-n);  // would be B-t for left tail
  }

}

//=====================================================================

Double_t myCrystalBallLeftTail(double* x, double* par) {

  Double_t xcur = x[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t invAbsAlpha = 1./absAlpha;

  if ( t >= -absAlpha)  {   // would be t <= absAlpha for right tail
    return N*exp(-0.5*t*t);
  } else {
    Double_t A = TMath::Power(n*invAbsAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t B = n*invAbsAlpha - absAlpha;
    return N*A * TMath::Power(B-t,-n);  // would be B+t for right tail
  }

}

//=====================================================================

Double_t my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Double_t xcur = x[0];
  Double_t alphaL = par[0];
  Double_t nL = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t alphaR = par[5];
  Double_t nR = par[6];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t invAbsAlphaL = 1./absAlphaL;
  Double_t absAlphaR = fabs((Double_t)alphaR);
  Double_t invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

//=====================================================================

void setHistColor(vector<Int_t> &histColor, const Int_t nObject) {
  
  Int_t colorList[] = {kBlue, kRed, kGreen, kOrange+1, kCyan, kViolet};  
  // the first color is for the main object. This array may contain more values than nSamples  
  Int_t vectorIsEmpty = ((Int_t)histColor.size()) == 0 ? 1 : 0;
                                                          
  for (Int_t i = 0; i < nObject; i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)         

    if (vectorIsEmpty) histColor.push_back(colorList[i]);
    else histColor[i] = colorList[i];    // if histColor is just uninitialized but not really empty

  }
  
}


void buildChainWithFriend(TChain* chain, TChain* chFriend, string sampleName) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA") {

    // 2016 2.6fb^-1
    if (DATA2016) {
      if (SKIM_1LEP1JET_80X) {
	subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v1_runs_272023_273146");
	subSampleNameVector.push_back("SingleElectron_Run2016B_PromptReco_v2_runs_273150_275376");
	subSampleNameVector.push_back("SingleElectron_Run2016C_PromptReco_v2_runs_275420_276097");
      } else {
	subSampleNameVector.push_back("SingleElectron_PromptReco_v1_runs_272021_273149");
	subSampleNameVector.push_back("SingleElectron_PromptReco_v2_runs_273150_274443");
      }
    } else {
    // 2015 2.32 fb^-1
      subSampleNameVector.push_back("SingleElectron_Run2015C_16Dec_runs_254227_254914");
      subSampleNameVector.push_back("SingleElectron_Run2015D_16Dec_runs_256630_260627");
    }

  } else if (sampleName == "WJetsToLNu") {

    if (SKIM_1LEP1JET_80X) {
      subSampleNameVector.push_back("WJetsToLNu_HT100to200");
      subSampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT200to400");
      subSampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT400to600");
      subSampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT600to800");
      subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
      subSampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
      subSampleNameVector.push_back("WJetsToLNu_HT1200to2500");
      subSampleNameVector.push_back("WJetsToLNu_HT2500toInf");
    } else {    
      subSampleNameVector.push_back("WJetsToLNu_HT100to200");
      subSampleNameVector.push_back("WJetsToLNu_HT200to400");
      subSampleNameVector.push_back("WJetsToLNu_HT400to600");
      subSampleNameVector.push_back("WJetsToLNu_HT600to800");
      subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
      if (!DATA2016) subSampleNameVector.push_back("WJetsToLNu_HT1200to2500"); // note, 1200to2500 bin missing for 2016 MC
      subSampleNameVector.push_back("WJetsToLNu_HT2500toInf"); 
    }
  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  //2016 trees
  string treePath = "";
  if (SKIM_1LEP1JET_80X) treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_1LEP1JET_80X/"; // 2016 trees with skim 1 lep1jet
  else {
    if (DATA2016) treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_1LEPSKIM_80X/"; // 2016 trees
    else treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_25ns_1LEPSKIM_76X/";   //2015 trees
  }

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    string friend_treeRootFile = "";
    if (SKIM_1LEP1JET_80X) {
      friend_treeRootFile = treePath + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    } else {
      if (DATA2016) friend_treeRootFile = treePath + "friends_evVarFriend_" + subSampleNameVector[i]+ ".root";
      else friend_treeRootFile = treePath + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    }

    chain->Add(TString(treeRootFile.c_str()));
    chFriend->Add(TString(friend_treeRootFile.c_str()));

  }

  cout << "Adding friend to chain ..." << endl;
  chain->AddFriend(chFriend);  //adding whole friend chain as friend                                                           

  if(!chain || !chFriend) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  cout << chain->GetEntries() << endl;

}

//============================================================   

Int_t getBinNumber(const Float_t value, const vector<Float_t> &binEdgesVector) {

  // return invalid value if bin not found
  // return -1 if value > binEdgesVector.last(),
  // return -2 if value < binEdgesVector.first(),

  if (value > binEdgesVector[0]) {

    Int_t bin = 0;
    Int_t howManyBins = binEdgesVector.size() - 1;

    while (bin < howManyBins) {
      if (value < binEdgesVector[bin+1]) {
	return bin;
      }
      bin++;
    }      

    return -1;

  } else return -2;

}

//============================================================

void EoverP::Loop(const string sampleName, const vector<Float_t> &corrEnergybinEdges, const string& dirName)
{


   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nEle40T",1);
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus                                               
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_correctedEcalEnergy",1);
   fChain->SetBranchStatus("LepGood_superCluster_rawEnergy",1);
   fChain->SetBranchStatus("LepGood_eSuperClusterOverP",1);
   fChain->SetBranchStatus("LepGood_r9",1);
   fChain->SetBranchStatus("met_pt",1);

   // branches for MC study
   if (sampleName != "DATA"){
     fChain->SetBranchStatus("ngenLep",1);
     fChain->SetBranchStatus("genLep_motherId",1);
     fChain->SetBranchStatus("genLep_pdgId",1);
     fChain->SetBranchStatus("genLep_pt",1);
     fChain->SetBranchStatus("genLep_eta",1);
     fChain->SetBranchStatus("genLep_phi",1);
     fChain->SetBranchStatus("genLep_mass",1);
   }

   gStyle->SetStatStyle(0);

   // passing from main
   // vector<Float_t> corrEnergybinEdges;
   // corrEnergybinEdges.push_back(25.0);
   // corrEnergybinEdges.push_back(50.0);
   // corrEnergybinEdges.push_back(75.0);
   // corrEnergybinEdges.push_back(100.0);
   // corrEnergybinEdges.push_back(150.0);
   // corrEnergybinEdges.push_back(200.0);
   // corrEnergybinEdges.push_back(250.0);
   // corrEnergybinEdges.push_back(300.0);
   // corrEnergybinEdges.push_back(350.0);
   // corrEnergybinEdges.push_back(450.0);
   // corrEnergybinEdges.push_back(600.0);
   // corrEnergybinEdges.push_back(750.0);
   // corrEnergybinEdges.push_back(900.0);

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

   Int_t nCorrEnergyBins = corrEnergybinEdges.size() -1;

   string rootfileName = dirName + "EoverP_" + sampleName + ".root";

   TFile *rootFile = new TFile((rootfileName).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout << "Error: file \"" << rootfileName << "\" was not opened." << endl;
     exit(EXIT_FAILURE);
   }

   vector<TH1F*> hEoverP_corrEnergyBin(nCorrEnergyBins,NULL);
   // histogram to store the mean energy in a given bin, useful to plot E/P asf of the mean energy in the bin with a TGraph
   TH1F* hMeanEnergyInCorrEnergyBin = new TH1F("hMeanEnergyInCorrEnergyBin","",nCorrEnergyBins,corrEnergybinEdges.data());

   vector<TH1F*> hEcorrOverEtrue_corrEnergyBin(nCorrEnergyBins,NULL);
   vector<TH1F*> hErawOverEtrue_corrEnergyBin(nCorrEnergyBins,NULL);
   vector<TH1F*> hPtrackOverEtrue_corrEnergyBin(nCorrEnergyBins,NULL);

   for (Int_t i = 0; i < nCorrEnergyBins; i++) {
     hEoverP_corrEnergyBin[i] = new TH1F(Form("hEoverP_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",200,0.05,2.05);
   }

   if (sampleName != "DATA") {
     for (Int_t i = 0; i < nCorrEnergyBins; i++) {
       hEcorrOverEtrue_corrEnergyBin[i] = new TH1F(Form("hEcorrOverEtrue_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",200,0.55,1.55);
       hErawOverEtrue_corrEnergyBin[i] = new TH1F(Form("hErawOverEtrue_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",200,0.55,1.55);
       hPtrackOverEtrue_corrEnergyBin[i] = new TH1F(Form("hPtrackOverEtrue_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",200,0.05,2.05);
     }
   }

   TH1F* hEoP_template = new TH1F("hEoP_template","template",200,0.05,2.05);

   // Float_t lowMeanCut = 0.85;
   // Float_t upMeanCut = 1.15;

   // used only for MC study
   /////////////////////////
   Int_t MCtruthMatchFound = 0;
   Float_t Etrue = 0.0; 
   /////////////////////////

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%500000 == 0) cout << jentry << endl;
  
      if (met_pt < 50) continue;
      if (!(nEle10V == 1 && nEle40T == 1)) continue;
      if (!(fabs(LepGood_pdgId[0]) == 11 && LepGood_pt[0] > 40 && fabs(LepGood_eta[0]) < 1.0 && LepGood_r9[0] > 0.94)) continue;
      //if (!(fabs(LepGood_pdgId[0]) == 11 && LepGood_pt[0] > 40 && fabs(LepGood_eta[0]) > 1.0 && fabs(LepGood_eta[0]) < 1.479)) continue;

      // here goes with the algorithm to match gen to reco electrons and in case fill histograms

      if (sampleName != "DATA") {

	MCtruthMatchFound = 0;  // reset to 0 
	Etrue = 0.0;            // reset this too
	Int_t i = 0;

	while (!MCtruthMatchFound && i < ngenLep) {

	  // start asking for e+/- (|pdgID| = 11) coming from a W+/- (|pdgID| = 24)
	    
	  if (fabs(genLep_pdgId[i]) == 11 && fabs(genLep_motherId[i]) == 24) {

	    TVector3 eleReco;
	    eleReco.SetPtEtaPhi(LepGood_pt[0],LepGood_eta[0],LepGood_phi[0]);
	    TVector3 eleGen;
	    eleGen.SetPtEtaPhi(genLep_pt[i],genLep_eta[i],genLep_phi[i]);

	    // now match is ok if deltaR < 0.3
	    // in this case, compute Etrue and set MCtruthMatchFound flag as 1

	    if (eleGen.DeltaR(eleReco) < 0.3) {

	      Etrue = eleGen.Mag();  // neglect mass since we have electrons
	      MCtruthMatchFound = 1;  

	    }
	    
	  }

	  i++;

	}

      } 

      // end of loop to match gen and reco electrons

      Double_t energyToUse = LepGood_correctedEcalEnergy[0];
      if (!USE_E) {
	Double_t theta = 2. * atan(exp(-LepGood_eta[0])); 
	energyToUse *= sin(theta);
      } 

      hEoP_template->Fill(LepGood_eSuperClusterOverP[0]);

      // look for the bin in the LepGood_correctedEcalEnergy variable      
      Int_t bin = getBinNumber(energyToUse,corrEnergybinEdges);  // this function returns negative value if bin not found

      if (bin >= 0) {

	hEoverP_corrEnergyBin[bin]->Fill(LepGood_eSuperClusterOverP[0]);

	//sum the energy to the bin content in the bin it belongs to (at the end we will divide by the number of entries in each bin)
	// using bin+1 because the histogram bin number goes from 1 to number of bins, while "bin" variable starts from 0
	hMeanEnergyInCorrEnergyBin->SetBinContent(bin+1, energyToUse + hMeanEnergyInCorrEnergyBin->GetBinContent(bin+1));  

	if (MCtruthMatchFound) {

	  hEcorrOverEtrue_corrEnergyBin[bin]->Fill(LepGood_correctedEcalEnergy[0]/Etrue);
	  hErawOverEtrue_corrEnergyBin[bin]->Fill(LepGood_superCluster_rawEnergy[0]/Etrue);
	  hPtrackOverEtrue_corrEnergyBin[bin]->Fill(LepGood_correctedEcalEnergy[0]/(LepGood_eSuperClusterOverP[0]*Etrue));    
	  // Ptrack/Etrue = corrE * (EoverP)^-1 / Etrue

	}

      }


   }  // end of loop on entries

   for (Int_t i = 0; i < nCorrEnergyBins; i++) {

     //hEoverP_corrEnergyBin[i]->SetStats(0);  // keep stat box in file, and remove in directly the function that plots it on canvas
     // keep the following axis settings
     hEoverP_corrEnergyBin[i]->GetXaxis()->SetTitle("E / P");
     hEoverP_corrEnergyBin[i]->GetXaxis()->SetTitleSize(0.06);
     hEoverP_corrEnergyBin[i]->GetXaxis()->SetTitleOffset(0.8);
     hEoverP_corrEnergyBin[i]->GetYaxis()->SetTitle("events");
     hEoverP_corrEnergyBin[i]->GetYaxis()->SetTitleSize(0.055);
     hEoverP_corrEnergyBin[i]->GetYaxis()->SetTitleOffset(0.8);

     hMeanEnergyInCorrEnergyBin->SetBinContent(i+1, hMeanEnergyInCorrEnergyBin->GetBinContent(i+1)/hEoverP_corrEnergyBin[i]->GetEntries());

   }

   if (USE_E) {
     hMeanEnergyInCorrEnergyBin->GetXaxis()->SetTitle("corrected E [GeV]");
     hMeanEnergyInCorrEnergyBin->GetYaxis()->SetTitle("mean E [GeV]");
   } else {
     hMeanEnergyInCorrEnergyBin->GetXaxis()->SetTitle("corrected E_{T} [GeV]");
     hMeanEnergyInCorrEnergyBin->GetYaxis()->SetTitle("mean E_{T} [GeV]");
   }
   hMeanEnergyInCorrEnergyBin->GetXaxis()->SetTitleSize(0.06);
   hMeanEnergyInCorrEnergyBin->GetXaxis()->SetTitleOffset(0.8);
   hMeanEnergyInCorrEnergyBin->GetYaxis()->SetTitleSize(0.055);
   hMeanEnergyInCorrEnergyBin->GetYaxis()->SetTitleOffset(0.8);

   if (sampleName != "DATA") {

     for (Int_t i = 0; i < nCorrEnergyBins; i++) {

       hEcorrOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitle("E_{corr} / E_{true}");
       hEcorrOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitleSize(0.06);
       hEcorrOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitleOffset(0.75);
       hEcorrOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitle("events");
       hEcorrOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitleSize(0.055);
       hEcorrOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitleOffset(0.8);

       hErawOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitle("E_{raw} / E_{true}");
       hErawOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitleSize(0.06);
       hErawOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitleOffset(0.75);
       hErawOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitle("events");
       hErawOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitleSize(0.055);
       hErawOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitleOffset(0.8);

       hPtrackOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitle("P_{track} / E_{true}");
       hPtrackOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitleSize(0.06);
       hPtrackOverEtrue_corrEnergyBin[i]->GetXaxis()->SetTitleOffset(0.75);
       hPtrackOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitle("events");
       hPtrackOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitleSize(0.055);
       hPtrackOverEtrue_corrEnergyBin[i]->GetYaxis()->SetTitleOffset(0.8);
    
     }

   }

   rootFile->Write();
   rootFile->Close();
   delete rootFile;


}

void getTexMCSampleName(const string &MCSampleName, string &texMCSampleName) {
  
  if (MCSampleName == "WJetsToLNu") texMCSampleName = "W(l#nu)+jets";
  else if (MCSampleName == "DYJetsToLL") texMCSampleName = "Z(ll)+jets";

}

void plotDistribution(const string &sampleName, const vector<Float_t> &corrEnergybinEdges, TH1F* hPeak, TH1F* hSigma, const string &hNameID, const string &dirName, TH1F* hPeakShift_MCdata = NULL) {

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 
  TVirtualFitter::SetDefaultFitter("Minuit");

  string fileName = dirName + "EoverP_" + sampleName + ".root";

  TCanvas *c = new TCanvas("c",""); 

  TH1F* htmp = NULL;
  TH1F* hist = NULL;    

  TFile* f = TFile::Open(fileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH1F* hEoP_template = NULL;
  // get E/P template from MC file

  // -----------------------------------------------------------
  // get MC template. If running on data, open MC file, else the file is already the right one and I just need to get histograms

  TFile* fMC = NULL;  // to get MC file when running on data
  string fileNameMC = dirName + "EoverP_WJetsToLNu.root";  // template is created when running on MC, and is stored in this file

  if (sampleName == "DATA") {

    //    getEoP_templateMC(hEoP_template, dirName);
    cout << "Running on data. Opening file " << fileNameMC << " to get template." << endl;

    fMC = TFile::Open(fileNameMC.c_str(),"READ");
    if (!fMC || !fMC->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<fileNameMC<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }
    
  //   // get E/P template from file
  //   htmp = (TH1F*)fMC->Get("hEoP_template");  
  //   if (!htmp) {
  //     cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileNameMC << "'. End of programme." << endl;
  //     exit(EXIT_FAILURE);
  //   }
  //   hEoP_template = (TH1F*) htmp->Clone();
    
  } else {

  //   htmp = (TH1F*)f->Get("hEoP_template");  
  //   if (!htmp) {
  //     cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileName << "'. End of programme." << endl;
  //     exit(EXIT_FAILURE);
  //   }
  //   hEoP_template = (TH1F*) htmp->Clone();

  // }

  // if (!hEoP_template) {
  //   cout << "Error: hEoP_template is NULL! End of programme." << endl;
  //   exit(EXIT_FAILURE);
  }

  // end of access of MC file to get templates
  // -----------------------------------------------------------

  UInt_t nBins = corrEnergybinEdges.size() -1;

  for (UInt_t i = 0; i < nBins; i++) {

    htmp = (TH1F*)f->Get(Form("h%s_corrEnergyBin%1.0fTo%1.0f",hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1]));  
    if (!htmp) {
      cout << "Error: histogram " << htmp->GetName() << " not found in file ' " << fileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    // now, if using data I will fit with template obtained from the same distribution in MC file. Since the histogram's name is the same in both file, I 
    // got a strange behaviour, that is, when using fMC->Get(<name>) I actually get the same histogram in data. A solution is deleting the previous histogram
    // got from data file before getting that for MC. In this case, better to give the clone a different name to make explicit it is for data 
    if (sampleName == "DATA") hist = (TH1F*)htmp->Clone(Form("%s_data",htmp->GetName()));
    else hist = (TH1F*)htmp->Clone();

    //    hist->SetStats(0);
    gStyle->SetOptStat(10);
    gStyle->SetOptFit(0112);
    hist->Draw("HE");
    if (sampleName == "DATA") {
      if (corrEnergybinEdges[i] > 449.9) hist->Rebin(5); 
      else if (corrEnergybinEdges[i] > 349.9) hist->Rebin(4);
      else if (corrEnergybinEdges[i] > 274.9) hist->Rebin(3); 
      else if (corrEnergybinEdges[i] > 174.9) hist->Rebin(2); 
      //else if (corrEnergybinEdges[i] > 199.9) hist->Rebin(2); // when using E instead of Et 249.9 is ok
    } else if (hNameID != "EoverP") {
      if (hNameID == "PtrackOverEtrue") {
	if (corrEnergybinEdges[i] > 449.9) hist->Rebin(4); 
	else if (corrEnergybinEdges[i] > 349.9) hist->Rebin(3); 
	else if (corrEnergybinEdges[i] > 174.9) hist->Rebin(2); 
      } else {
	hist->GetXaxis()->SetRangeUser(0.85,1.15);
      }
    } else {
      if (corrEnergybinEdges[i] > 449.9) hist->Rebin(3); 
      else if (corrEnergybinEdges[i] > 224.9) hist->Rebin(2);      
    }
    c->Update();

    // fitting:
    // do a first fit with a simple gaussian in the core
    Double_t gaussEdgeL = 0.9;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
    Double_t gaussEdgeR = 1.1;  //right side ...
    if (corrEnergybinEdges[i] > 349.9) {
      gaussEdgeL = 0.75;
      gaussEdgeR = 1.25;
    }
    if (sampleName != "DATA" && hNameID != "EoverP") {
      if (hNameID == "PtrackOverEtrue") {
	gaussEdgeL = 0.9;
	gaussEdgeR = 1.1;
      } else {
	gaussEdgeL = 0.95;
	gaussEdgeR = 1.05;	
      }
    }
    hist->Fit("gaus","E L I Q 0","",gaussEdgeL,gaussEdgeR);  // L: loglikelihood method, 0: do not plot this fit, Q: quiet mode (minimum printing)
    Double_t gaussNorm = hist->GetFunction("gaus")->GetParameter(0);
    Double_t gaussMean = hist->GetFunction("gaus")->GetParameter(1);
    //Double_t gaussMeanError = hist->GetFunction("gaus")->GetParError(1);
    Double_t gaussSigma = hist->GetFunction("gaus")->GetParameter(2);
    // now use crystal ball with right tail
    Double_t funcRangeLeft = gaussEdgeL;
    Double_t funcRangeRight = 2.0;
    if (FIT_2SIDE_CB) {
      if (hNameID != "EoverP") {
	if (hNameID != "PtrackOverEtrue") {
	  funcRangeLeft = 0.8;
	  funcRangeRight = gaussEdgeR;
	} else {
	  funcRangeLeft = 0.2;
	  funcRangeRight = 1.3;
	}
      } else {
	funcRangeLeft = 0.5;
      }
    }
    else {
      if (hNameID != "EoverP") {
	if (hNameID != "PtrackOverEtrue") funcRangeLeft = 0.8;
	else funcRangeLeft = 0.2;
	funcRangeRight = gaussEdgeR;
      }
    }

    TF1 *cb1;

    if (FIT_2SIDE_CB) {

      cb1 = new TF1("cb1",&my2sideCrystalBall,funcRangeLeft,funcRangeRight,7);
      cb1->SetParNames("alphaL","nL","mu","sigma","N","alphaR","nR");  
      cb1->SetParLimits(cb1->GetParNumber("nL"),0.1,15);
      cb1->SetParLimits(cb1->GetParNumber("nR"),0.1,15);
      cb1->SetParLimits(cb1->GetParNumber("alphaL"),-10.,-0.01); 
      cb1->SetParLimits(cb1->GetParNumber("alphaR"),0.01,10);
      cb1->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm,(gaussEdgeR-gaussMean)/gaussSigma,5);
      cb1->SetLineColor(kRed);

    } else {

      if (hNameID != "EoverP") cb1 = new TF1("cb1",&myCrystalBallLeftTail,funcRangeLeft,funcRangeRight,5);  // last parameter is the number of free parameters
      else cb1 = new TF1("cb1",&myCrystalBallRightTail,funcRangeLeft,funcRangeRight,5);
      cb1->SetParNames("alpha","n","mu","sigma","N");  
      cb1->SetParLimits(cb1->GetParNumber("n"),0.1,15); 
      if (hNameID != "EoverP") {
	cb1->SetParLimits(cb1->GetParNumber("alpha"),-10.0,-0.01);
	cb1->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm);
      } else {
	cb1->SetParLimits(cb1->GetParNumber("alpha"),0.01,10);
	cb1->SetParameters((gaussEdgeR-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm);
      }
      // with the following (or some of them) it looks like the fit doesn't work well, tipically the value from fit is out of the range
      // cb1->SetParLimits(cb1->GetParNumber("mu"),gaussMean-3.0*gaussSigma,gaussMean+3.0*gaussSigma);
      // cb1->SetParLimits(cb1->GetParNumber("sigma"),0.1*gaussSigma,10.0*gaussSigma);
      // cb1->SetParLimits(cb1->GetParNumber("N"),0.1*gaussNorm,10.0*gaussNorm);
    }

    // build TLegend
    TLegend *leg = new TLegend(0.11,0.8,0.4,0.89);
    TLegend *legFitFunction = new TLegend(0.11,0.5,0.35,0.8);

    // do the fit
    // before fitting, clone hist . Will use the clone to do the fit with the template. This is needed to plot 2 statistic boxes for crystal Ball and template fit
    TH1F* histClone = (TH1F*) hist->Clone("histClone");
    TFitResultPtr frp1 = hist->Fit(cb1,"E L I S Q B R","HE",funcRangeLeft,funcRangeRight);
    leg->AddEntry(hist,"distribution","l");
    legFitFunction->SetHeader("fit functions:");
    legFitFunction->AddEntry(cb1,"Crystal Ball","l");
    TFitResultPtr frp_template;
    TF1* f_EoP_template = NULL;

    if (sampleName == "DATA" && hNameID == "EoverP") {

      // ----------------------------------------------------------------
      // fitting with template (only for data and for E/P
      
      // ------------------------------------------------------------
      // get binned MC histograms to use them as template
	
      //    getEoP_templateMC(hEoP_template, dirName);
    
      // get E/P template from file

      // before getting the template, delete the same named histogram from data file, otherwise I am getting the same data histogram despite the fact that
      // I'm using the pointer to MC file
      delete htmp;
      fMC->cd();
      TH1F* htmpMC = NULL;
      htmpMC = (TH1F*)fMC->Get(Form("hEoverP_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));  
      if (!htmpMC) {
	cout << "Error: histogram " << htmpMC->GetName() << " not found in file ' " << fileNameMC << "'. End of programme." << endl;
	exit(EXIT_FAILURE);
      } else {
	//cout << "Energy ["<< corrEnergybinEdges[i] <<"," <<corrEnergybinEdges[i+1]<< "] \t max htmpMC: " << htmpMC->GetBinContent(htmpMC->GetMaximumBin()) << endl;
      }

      hEoP_template = (TH1F*) htmpMC->Clone(Form("hEoverP_corrEnergyBin%1.0fTo%1.0f_MC",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));
      if (!hEoP_template) {
	cout << "Error: hEoP_template is NULL! End of programme." << endl;
	exit(EXIT_FAILURE);
      }	
      // rebinning if needed
      if (corrEnergybinEdges[i] > 449.9) hEoP_template->Rebin(3); 
      else if (corrEnergybinEdges[i] > 174.9) hEoP_template->Rebin(2);      

      // end of access of MC file to get templates
      // -----------------------------------------------------------

      // class defining the template function
      histoFunc* templateHistoFunc = new histoFunc(hEoP_template); 
      //cout << "templateHistoFunc->GetIntegral() : " << templateHistoFunc->GetIntegral() << endl;
      Double_t templateFitRangeLow = 0.6;  // 0.7
      Double_t templateFitRangeUp = 2.0;   // 2.0 
      f_EoP_template = new TF1("f_EoP_template", templateHistoFunc, templateFitRangeLow, templateFitRangeUp, 3, "histoFunc");
      f_EoP_template -> SetParName(0,"Norm"); 
      f_EoP_template -> SetParName(1,"Scale factor");
      f_EoP_template -> SetParName(2,"Shift of x"); 
      //      f_EoP_template -> SetLineWidth(1); 
      f_EoP_template -> SetNpx(10000);
      //Double_t xNorm = histClone->GetEntries()/hEoP_template->GetEntries() * histClone->GetBinWidth(1)/hEoP_template->GetBinWidth(1); // GetEntries includes under/overflow
      Double_t xNorm = histClone->Integral(getLowBinGivenRange(histClone,templateFitRangeLow),getUpBinGivenRange(histClone,templateFitRangeUp)) /
	hEoP_template->Integral(getLowBinGivenRange(hEoP_template,templateFitRangeLow),getUpBinGivenRange(hEoP_template,templateFitRangeUp)) * 
	histClone->GetBinWidth(1)/hEoP_template->GetBinWidth(1); 
      // f_EoP_template -> SetParameter(0, xNorm);
      // f_EoP_template -> SetParLimits(0, 0.95 * xNorm, 1.05 * xNorm);
      f_EoP_template -> FixParameter(0, xNorm);
      f_EoP_template -> SetParameter(1, 1.0);  // maybe this should be the inverse of the width ratio
      f_EoP_template -> SetParLimits(1,0.1,10.0); 
      Double_t MC_data_peakDeltaX = hEoP_template->GetBinCenter(hEoP_template->GetMaximumBin()) - frp1->Parameter(2);
      f_EoP_template -> SetParameter(2, MC_data_peakDeltaX);  // will set this to Dx between peak in data and MC ( MC - data )
      f_EoP_template -> SetParLimits(2, 0.1 * MC_data_peakDeltaX, 10.0 * MC_data_peakDeltaX);  // will set this to Dx between peak in data and MC ( MC - data )
      //Double_t MC_data_peakDeltaX = 0.2;
      //f_EoP_template -> FixParameter(2, MC_data_peakDeltaX);  // will set this to Dx between peak in data and MC ( MC - data )
      f_EoP_template -> SetLineColor(kGreen+2);
      // -----------------------------------------------------------------

      frp_template = histClone->Fit("f_EoP_template","E L I S Q B R","HE SAMES", templateFitRangeLow,templateFitRangeUp);
      legFitFunction->AddEntry(f_EoP_template,"MC template","l");

    }

    leg->Draw();
    leg->SetMargin(0.3); 
    leg->SetBorderSize(0);

    legFitFunction->Draw();
    legFitFunction->SetMargin(0.3); 
    legFitFunction->SetBorderSize(0);

    c->Update();
    // box for fit with Crystal Ball
    TPaveStats *stat = (TPaveStats*)(hist->FindObject("stats"));
    if(stat) {
      // stat->SetTextColor(kBlue);
      // stat1->SetTextColor(kGreen);
      float width = stat->GetX2NDC() - stat->GetX1NDC();
      // make stat box bigger
      stat->SetX1NDC(stat->GetX1NDC() - 0.5 * width);     
      stat->Draw();
    }
    // Box for fit with template (data only)
    TPaveStats *stat1 = NULL;
    if (sampleName == "DATA") {
      histClone->Draw("HE SAME");
      c->Update();
      stat1 = (TPaveStats*)(histClone->FindObject("stats"));
      if (stat1) {
	float height = stat1->GetY2NDC() - stat1->GetY1NDC();
	stat1->SetY1NDC(stat->GetY1NDC() - height);
	stat1->SetY2NDC(stat->GetY1NDC() );
	stat1->SetX1NDC(stat->GetX1NDC());
	stat1->Draw();
      }
    }

    if (hPeak != NULL) {
      hPeak->SetBinContent(i+1, frp1->Parameter(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
      hPeak->SetBinError(i+1, frp1->ParError(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
    }
    if (hPeakShift_MCdata  != NULL) {
      hPeakShift_MCdata ->SetBinContent(i+1, frp_template->Parameter(2)); // 2 is deltaX 
      hPeakShift_MCdata ->SetBinError(i+1, frp_template->ParError(2)); // 2 is deltaX 
    }
    if (hSigma != NULL) {
      hSigma->SetBinContent(i+1, frp1->Parameter(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
      hSigma->SetBinError(i+1, frp1->ParError(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
    }

    if (USE_E) {      
      hist->SetTitle(Form("%1.0f < E [GeV] < %1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_E%1.0fTo%1.0f_%s.pdf",dirName.c_str(),hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_E%1.0fTo%1.0f_%s.png",dirName.c_str(),hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));
    } else {
      hist->SetTitle(Form("%1.0f < E_{T}[GeV] < %1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));
      c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.pdf",dirName.c_str(),hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));
      c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.png",dirName.c_str(),hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));
    }

  }

  //cout << "CHECK 1" << endl;
  delete c;
  //cout << "CHECK 2" << endl;
  //delete htmp;
  //cout << "CHECK 3" << endl;
  //delete hist;  // if I delete hist I get segmentation fault (it didn't happen before trying to fit with template)
  //cout << "CHECK 4" << endl;

}

//==========================================================================================


void drawPlotDataMC(TH1F* hdata, TH1F* hmc, const string& MCSampleName, const string& xAxisName, const string& yAxisName, const string& canvasName, const string &dirName) {

  TCanvas *cfit = new TCanvas("cfit","",700,700);

  // now here we go with the canvas                                                                                                                                    
  TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas                                                                                           
  TPad *subpad_2 = NULL;
  TLegend *leg = new TLegend(0.5,0.7,0.79,0.90);

  subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  subpad_2->SetGridy();
  subpad_2->SetBottomMargin(0.3);
  subpad_1->Draw();
  subpad_2->Draw();
  subpad_1->cd();
  
  hdata->SetStats(0);
  hdata->SetLineColor(kRed);
  hdata->Draw("HE");
  hdata->GetXaxis()->SetLabelSize(0.45);
  hdata->GetYaxis()->SetTitle(yAxisName.c_str());
  hdata->GetYaxis()->SetTitleSize(0.06);
  hdata->GetYaxis()->SetTitleOffset(0.8);
  if (SET_SCALE_ON_Y) {
    if (yAxisName == "peak(E/P)") {
      hdata->SetMinimum(0.94);
      hdata->SetMinimum(0.965);
      hdata->SetMaximum(1.005);
    } else if (yAxisName == "#sigma(E/P)") {
      hdata->SetMinimum(0.02);
      hdata->SetMaximum(0.15);
    }
  }
  hmc->SetLineColor(kBlue);
  hmc->Draw("HE SAME");

  string texMCSampleName = "";
  getTexMCSampleName(MCSampleName, texMCSampleName);

  leg->AddEntry(hdata,"data","lf");
  leg->AddEntry(hmc,Form("%s",texMCSampleName.c_str()),"lf");
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend

  //ratio plot

  subpad_2->cd();
  TH1F * ratioplot = NULL; // will use it for the ratio plots                                                                                                   
  ratioplot = new TH1F(*hdata);  // to have ratioplot with the same x axis range as hdata (or hmc), I copy it from hdata created above, then I substitute bin content with hdata/hmc                                                                                                                                                   
  ratioplot->Divide(hdata,hmc);
  ratioplot->SetStats(0);
  ratioplot->SetTitle("");
  ratioplot->GetXaxis()->SetLabelSize(0.10);
  ratioplot->GetXaxis()->SetTitle(xAxisName.c_str());
  ratioplot->GetXaxis()->SetTitleSize(0.14);
  ratioplot->GetXaxis()->SetTitleOffset(0.8);
  ratioplot->GetYaxis()->SetLabelSize(0.10);
  ratioplot->GetYaxis()->SetTitle("data / MC");
  ratioplot->GetYaxis()->SetTitleSize(0.12);
  ratioplot->GetYaxis()->SetTitleOffset(0.45);
  ratioplot->GetYaxis()->CenterTitle();
  ratioplot->GetYaxis()->SetNdivisions(011);
  ratioplot->SetMarkerStyle(8);  //medium dot  
  TH1F *ratioplotCopy = (TH1F*) ratioplot->DrawCopy("E");
  //ratioplot->Draw("E");
  if (SET_SCALE_ON_Y) {
    if (yAxisName == "peak(E/P)") {
      ratioplotCopy->SetMinimum(0.97);
      ratioplotCopy->SetMinimum(0.985);
      ratioplotCopy->SetMaximum(1.015);
    } else if (yAxisName == "#sigma(E/P)") {
      // ratioplotCopy->SetMinimum(0.9);
      // ratioplotCopy->SetMaximum(0.14);
    }
  }

  cfit->SaveAs( (dirName + canvasName + ".pdf").c_str() );
  cfit->SaveAs( (dirName + canvasName + ".png").c_str() );

  delete cfit;

}

//==========================================================================================


void drawPlotOnlyMC(vector<TH1F*> &hmcVector, const vector<string> &legEntryName, const string& MCSampleName, const string& xAxisName, const string& yAxisName, const string& canvasName, const string &dirName) {

  TCanvas *cfitMC = new TCanvas("cfitMC","");

  // now here we go with the canvas                                                                                                                                    
  // TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas                                                                                        
  // TPad *subpad_2 = NULL;
  TLegend *leg = new TLegend(0.5,0.6,0.79,0.90);

  // subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
  // subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
  // subpad_2->SetGridy();
  // subpad_2->SetBottomMargin(0.3);
  // subpad_1->Draw();
  // subpad_2->Draw();
  // subpad_1->cd();

  vector<Int_t> histColor;
  setHistColor(histColor,(Int_t)hmcVector.size());
  
  Double_t maximumYaxisValue = 0.0;
  Double_t minimumYaxisValue = 100000.0;
  
  for (UInt_t i = 0; i < hmcVector.size(); i++) {

    Double_t value = hmcVector[i]->GetBinContent(hmcVector[i]->GetMaximumBin()) + hmcVector[i]->GetBinError(hmcVector[i]->GetMaximumBin()); 
    if ( value > maximumYaxisValue) maximumYaxisValue = value; 
    value = hmcVector[i]->GetBinContent(hmcVector[i]->GetMinimumBin()) - hmcVector[i]->GetBinError(hmcVector[i]->GetMinimumBin());
    if ( value < minimumYaxisValue) minimumYaxisValue = value;


  } 

  for (UInt_t i = 0; i < hmcVector.size(); i++) {

    hmcVector[i]->SetStats(0);
    hmcVector[i]->SetLineColor(histColor[i]);
    if ( i == 0 ) {
      hmcVector[i]->SetMaximum(maximumYaxisValue * 1.02);  // slightly increase the y scale for maximum
      hmcVector[i]->SetMinimum(minimumYaxisValue * 0.98);  // slightly decrease the y scale for minimum
      if (hmcVector[i]->GetMinimum() < 0.02) hmcVector[i]->SetMinimum(0.0);
      hmcVector[i]->Draw("HE");
      hmcVector[i]->GetXaxis()->SetTitle(xAxisName.c_str());
      hmcVector[i]->GetXaxis()->SetLabelSize(0.05);
      hmcVector[i]->GetXaxis()->SetTitleSize(0.06);
      hmcVector[i]->GetXaxis()->SetTitleOffset(0.8);
      hmcVector[i]->GetYaxis()->SetTitle(yAxisName.c_str());
      hmcVector[i]->GetYaxis()->SetTitleSize(0.06);
      hmcVector[i]->GetYaxis()->SetTitleOffset(0.8);
      if (SET_SCALE_ON_Y) {
	if (yAxisName == "peak position") hmcVector[i]->GetYaxis()->SetRangeUser(0.98,1.04);
	else if (yAxisName == "#sigma of distribution") hmcVector[i]->GetYaxis()->SetRangeUser(0.0,0.11);
      }
    } else {
      hmcVector[i]->Draw("HE SAME");
    }

  }

  string texMCSampleName = "";
  getTexMCSampleName(MCSampleName, texMCSampleName);

  leg->AddEntry((TObject*)0,texMCSampleName.c_str(),"");
  for (UInt_t i = 0; i < hmcVector.size(); i++) {
    leg->AddEntry(hmcVector[i],Form("%s",legEntryName[i].c_str()),"lf");
  }
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend

  //ratio plot

  // subpad_2->cd();
  // TH1F * ratioplot = NULL; // will use it for the ratio plots                                                                                                   
  // ratioplot = new TH1F(*hdata);  // to have ratioplot with the same x axis range as hdata (or hmc), I copy it from hdata created above, then I substitute bin content with hdata/hmc                                                                                                                                                   
  // ratioplot->Divide(hdata,hmc);
  // ratioplot->SetStats(0);
  // ratioplot->SetTitle("");
  // ratioplot->GetXaxis()->SetLabelSize(0.10);
  // ratioplot->GetXaxis()->SetTitle(xAxisName.c_str());
  // ratioplot->GetXaxis()->SetTitleSize(0.14);
  // ratioplot->GetXaxis()->SetTitleOffset(0.8);
  // ratioplot->GetYaxis()->SetLabelSize(0.10);
  // ratioplot->GetYaxis()->SetTitle("data / MC");
  // ratioplot->GetYaxis()->SetTitleSize(0.12);
  // ratioplot->GetYaxis()->SetTitleOffset(0.45);
  // ratioplot->GetYaxis()->CenterTitle();
  // ratioplot->GetYaxis()->SetNdivisions(011);
  // ratioplot->SetMarkerStyle(8);  //medium dot  
  // ratioplot->DrawCopy("E");
  cfitMC->SaveAs( (dirName + canvasName + ".pdf").c_str() );
  cfitMC->SaveAs( (dirName + canvasName + ".png").c_str() );

  delete cfitMC;

}

//=================================================================================

void plotFromFit(const string &dataSampleName, const string &MCSampleName, const vector<Float_t> &corrEnergybinEdges, const string & dirName) {
  
  Int_t nCorrEnergyBins = corrEnergybinEdges.size() -1;

  TH1F* hPeakEoverPdata = new TH1F("hPeakEoverPdata","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEoverPdata = new TH1F("hSigmaEoverPdata","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hPeakEoverPmc = new TH1F("hPeakEoverPmc","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEoverPmc = new TH1F("hSigmaEoverPmc","",nCorrEnergyBins,corrEnergybinEdges.data());

  TH1F* hPeakShift_MCdata = new TH1F("hPeakShift_MCdata","",nCorrEnergyBins,corrEnergybinEdges.data()); 

  plotDistribution(dataSampleName, corrEnergybinEdges, hPeakEoverPdata, hSigmaEoverPdata, "EoverP",dirName, hPeakShift_MCdata);
  plotDistribution(MCSampleName, corrEnergybinEdges, hPeakEoverPmc, hSigmaEoverPmc, "EoverP",dirName);

  TH1F* hPeakEoverPdata_templateFit = NULL;
  if (hPeakShift_MCdata) hPeakEoverPdata_templateFit = new TH1F("hPeakEoverPdata_templateFit","",nCorrEnergyBins,corrEnergybinEdges.data());
  if (hPeakEoverPdata_templateFit) {
    hPeakEoverPdata_templateFit->Sumw2();
    hPeakEoverPdata_templateFit->Add(hPeakEoverPmc,hPeakShift_MCdata,1.0,-1.0);  // data = MC - ("MC-Data")
  }

  if (USE_E) {
    drawPlotDataMC(hPeakEoverPdata, hPeakEoverPmc, MCSampleName, "corrected E [GeV]", "peak(E/P)", "modeEoverPfromFit",dirName);
    if (hPeakEoverPdata_templateFit) drawPlotDataMC(hPeakEoverPdata_templateFit, hPeakEoverPmc, MCSampleName, "corrected E [GeV]", "peak(E/P)", "modeEoverPfromTemplateFit",dirName);
    drawPlotDataMC(hSigmaEoverPdata, hSigmaEoverPmc, MCSampleName, "corrected E [GeV]", "#sigma(E/P)", "sigmaEoverPfromFit",dirName);
  } else {
    drawPlotDataMC(hPeakEoverPdata, hPeakEoverPmc, MCSampleName, "corrected E_{T} [GeV]", "peak(E/P)", "modeEoverPfromFit",dirName);
    if (hPeakEoverPdata_templateFit) drawPlotDataMC(hPeakEoverPdata_templateFit, hPeakEoverPmc, MCSampleName, "corrected E_{T} [GeV]", "peak(E/P)", "modeEoverPfromTemplateFit",dirName);
    drawPlotDataMC(hSigmaEoverPdata, hSigmaEoverPmc, MCSampleName, "corrected E_{T} [GeV]", "#sigma(E/P)", "sigmaEoverPfromFit",dirName);
  }

  // quick plot of the shift between data and MC peak, as obtained from fit of data distribution with MC template

  TCanvas *cshift = new TCanvas("cshift","");
  if (hPeakShift_MCdata) {
    hPeakShift_MCdata->SetStats(0);
    hPeakShift_MCdata->Draw("HE");
    if (USE_E) hPeakShift_MCdata->GetXaxis()->SetTitle("corrected E [GeV]");
    else hPeakShift_MCdata->GetXaxis()->SetTitle("corrected E_{T} [GeV]");
    hPeakShift_MCdata->GetYaxis()->SetTitle("peak shift (MC - DATA) from template fit");
    hPeakShift_MCdata->GetXaxis()->SetLabelSize(0.05);
    hPeakShift_MCdata->GetXaxis()->SetTitleSize(0.06);
    hPeakShift_MCdata->GetXaxis()->SetTitleOffset(0.8);
    hPeakShift_MCdata->GetYaxis()->SetTitleSize(0.045);
    hPeakShift_MCdata->GetYaxis()->SetTitleOffset(1.1);
    cshift->SaveAs((dirName + "EoverP_peakShift_MCdata.png").c_str());
    cshift->SaveAs((dirName + "EoverP_peakShift_MCdata.pdf").c_str());
  }
  delete cshift;  

  // MC only study

  TH1F* hPeakEcorrOverEtrue = new TH1F("hPeakEcorrOverEtrue","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEcorrOverEtrue = new TH1F("hSigmaEcorrOverEtrue","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hPeakErawOverEtrue = new TH1F("hPeakErawOverEtrue","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaErawOverEtrue = new TH1F("hSigmaErawOverEtrue","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hPeakPtrackOverEtrue = new TH1F("hPeakPtrackOverEtrue","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaPtrackOverEtrue = new TH1F("hSigmaPtrackOverEtrue","",nCorrEnergyBins,corrEnergybinEdges.data());

  plotDistribution(MCSampleName, corrEnergybinEdges, hPeakEcorrOverEtrue, hSigmaEcorrOverEtrue, "EcorrOverEtrue",dirName);
  plotDistribution(MCSampleName, corrEnergybinEdges, hPeakErawOverEtrue, hSigmaErawOverEtrue, "ErawOverEtrue",dirName);
  plotDistribution(MCSampleName, corrEnergybinEdges, hPeakPtrackOverEtrue, hSigmaPtrackOverEtrue, "PtrackOverEtrue",dirName);    

  vector<TH1F*> hPeakVectorMC;
  hPeakVectorMC.push_back(hPeakEcorrOverEtrue);
  hPeakVectorMC.push_back(hPeakErawOverEtrue);
  hPeakVectorMC.push_back(hPeakPtrackOverEtrue);

  vector<TH1F*> hSigmaVectorMC;
  hSigmaVectorMC.push_back(hSigmaEcorrOverEtrue);
  hSigmaVectorMC.push_back(hSigmaErawOverEtrue);
  hSigmaVectorMC.push_back(hSigmaPtrackOverEtrue);

  vector<string> legEntryName;
  legEntryName.push_back("E_{corr}/E_{true}");
  legEntryName.push_back("E_{raw}/E_{true}");
  legEntryName.push_back("P_{track}/E_{true}");

  if (USE_E) {
    drawPlotOnlyMC(hPeakVectorMC, legEntryName, MCSampleName, "corrected E [GeV]", "peak position", "modeMCstudy",dirName); 
    drawPlotOnlyMC(hSigmaVectorMC, legEntryName, MCSampleName, "corrected E [GeV]", "#sigma of distribution", "sigmaMCstudy",dirName); 
  } else {
    drawPlotOnlyMC(hPeakVectorMC, legEntryName, MCSampleName, "corrected E_{T} [GeV]", "peak position", "modeMCstudy",dirName); 
    drawPlotOnlyMC(hSigmaVectorMC, legEntryName, MCSampleName, "corrected E_{T} [GeV]", "#sigma of distribution", "sigmaMCstudy",dirName); 
  }


}

//======================================================================


Int_t main(Int_t argc, char* argv[]) {

  Int_t doAll_flag = 1;
  Int_t doLoop_flag = 1;
  Int_t doPlot_flag = 1;

  string dirName = "";   // will be a path like "plot/<name>/" . Note the ending "/" 

  Option* option = new Option();
  option->Set_data2016(1);
  option->Set_skim1lep1jet80X(1);
  option->Set_fit2sideCB(0);
  option->Set_setScaleOnY(1);
  option->Set_useE(1);
 
  if (argc > 1) {

    for (Int_t i = 1; i < argc; i++) {

      string thisArgument(argv[i]);
      if (thisArgument == "-nl") {
	cout << "Passing option -nl: skip Loop on ntuples" << endl;
	doLoop_flag = 0;  // -nl --> no Loop
	doAll_flag = 0;
      } else if (thisArgument == "-np") {
	cout << "Passing option -np: skip creation of plots" << endl;
	doPlot_flag = 0;  // -np --> no plots
	doAll_flag = 0;
      } else if (thisArgument == "-dn") {   // -dn --> directory name
	cout << "Passing option -dn: passing name for directories to be created" << endl;
	dirName = string(argv[i+1]);
	option->Set_dirName(string(argv[i+1]));
	cout << "Saving output in '" << dirName << "'" << endl;
	i++;
      }

    }

  }

  vector<string> sampleName;
  sampleName.push_back("DATA");
  sampleName.push_back("WJetsToLNu");

  vector<Float_t> corrEnergybinEdges;
  // this binning looks ok
  corrEnergybinEdges.push_back(25.0);
  corrEnergybinEdges.push_back(50.0);
  corrEnergybinEdges.push_back(75.0);
  corrEnergybinEdges.push_back(100.0);
  corrEnergybinEdges.push_back(125.0);
  corrEnergybinEdges.push_back(175.0);
  corrEnergybinEdges.push_back(225.0);
  if (USE_E) {
    corrEnergybinEdges.push_back(275.0);
    corrEnergybinEdges.push_back(450.0);
    corrEnergybinEdges.push_back(900.0);
  } else {
    corrEnergybinEdges.push_back(275.0);
    //corrEnergybinEdges.push_back(350.0);
    corrEnergybinEdges.push_back(450.0);
    corrEnergybinEdges.push_back(900.0);
  }    

  // corrEnergybinEdges.push_back(25.0);
  // corrEnergybinEdges.push_back(50.0);
  // corrEnergybinEdges.push_back(75.0);
  // corrEnergybinEdges.push_back(100.0);
  // corrEnergybinEdges.push_back(125.0);
  // corrEnergybinEdges.push_back(150.0);
  // corrEnergybinEdges.push_back(200.0);
  // corrEnergybinEdges.push_back(250.0);
  // corrEnergybinEdges.push_back(325.0);
  // corrEnergybinEdges.push_back(450.0);
  // corrEnergybinEdges.push_back(900.0);


  if(doAll_flag || doLoop_flag) {

    for (UInt_t i = 0; i < sampleName.size(); i++) {
    
      // create chain                                                                                                         
    
      TChain* chain = new TChain("tree");
      TChain* chFriend = new TChain("mjvars/t");

      buildChainWithFriend(chain, chFriend, sampleName[i]);
    
      if(!chain || !chFriend) {
	cout << "Error: chain not created. End of programme" << endl;
	exit(EXIT_FAILURE);
      }

      EoverP tree(chain);
      tree.Loop(sampleName[i], corrEnergybinEdges, dirName);

      delete chain;
      delete chFriend;

    }

  }


  if(doAll_flag || doPlot_flag) {

    plotFromFit(sampleName[0],sampleName[1], corrEnergybinEdges, dirName);

  }

  // plotEoverP(sampleName[0],sampleName[1],"hMeanEoverP","< E/P >");
  // plotEoverP(sampleName[0],sampleName[1],"hModeEoverP"," mode of E/P");
  // plotEoverPdistribution(sampleName[0], corrEnergybinEdges);
  // plotEoverPdistribution(sampleName[1], corrEnergybinEdges);


  return 0;

}



////////////////////////////////////////
//TO BE IMPLEMENTED
//
//Plots of E/P asaf mean energy
//
