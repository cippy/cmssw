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

#define DATA2016 1

using namespace std;

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

// Double_t my2sideCrystalBall(string whichSide = "L",  double* x, double* par) {

//   //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

//   Double_t xcur = x[0];
//   Double_t alphaL = par[0];
//   Double_t nL = par[1];
//   Double_t mu = par[2];
//   Double_t sigma = par[3];
//   Double_t N = par[4];
//   Double_t alphaR = par[5];
//   Double_t nR = par[6];
//   Double_t t = (xcur-mu)/sigma;
//   Double_t absAlphaL = fabs((Double_t)alphaL);
//   Double_t invAbsAlphaL = 1./absAlphaL;
//   Double_t absAlphaR = fabs((Double_t)alphaR);
//   Double_t invAbsAlphaR = 1./absAlphaR;

  
//   if ( t<-absAlphaL ) {
//     //cout<<"checkpoint dscb left"<<endl;
//     Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
//     Double_t BL = nL*invAbsAlphaL - absAlphaL;
//     return N*AL*TMath::Power(BL-t,-nL);
//   } else if ( t <= absAlphaR )  {
//     //cout<<"checkpoint dscb gaussian"<<endl;
//     return N*exp(-0.5*t*t);
//   } else {
//     //cout<<"checkpoint dscb right"<<endl;
//     Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
//     Double_t BR = nR*invAbsAlphaR - absAlphaR;
//     return N*AR*TMath::Power(BR+t,-nR);
//   }

// }

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
      subSampleNameVector.push_back("SingleElectron_PromptReco_v1_runs_272021_273149");
      subSampleNameVector.push_back("SingleElectron_PromptReco_v2_runs_273150_274443");
    } else {
    // 2015 2.32 fb^-1
      subSampleNameVector.push_back("SingleElectron_Run2015C_16Dec_runs_254227_254914");
      subSampleNameVector.push_back("SingleElectron_Run2015D_16Dec_runs_256630_260627");
    }

  } else if (sampleName == "WJetsToLNu") {
    
    subSampleNameVector.push_back("WJetsToLNu_HT100to200");
    subSampleNameVector.push_back("WJetsToLNu_HT200to400");
    subSampleNameVector.push_back("WJetsToLNu_HT400to600");
    subSampleNameVector.push_back("WJetsToLNu_HT600to800");
    subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
    if (!DATA2016) subSampleNameVector.push_back("WJetsToLNu_HT1200to2500"); // note, 1200to2500 bin missing for 2016 MC
    subSampleNameVector.push_back("WJetsToLNu_HT2500toInf"); 

  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  //2016 trees
  string treePath = "";
  if (DATA2016) treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_1LEPSKIM_80X/"; // 2016 trees
  else treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_25ns_1LEPSKIM_76X/";   //2015 trees
  
  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    string friend_treeRootFile = "";
    if (DATA2016) friend_treeRootFile = treePath + "friends_evVarFriend_" + subSampleNameVector[i]+ ".root";
    else friend_treeRootFile = treePath + "evVarFriend_" + subSampleNameVector[i]+ ".root";

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

Double_t getEffectiveSigma(const TH1F* histo_) {
 
  //copied from Emanuele
 
  const TAxis *xaxis = histo_->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  // Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = histo_->GetMean();
  Double_t rms = histo_->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=histo_->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=histo_->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=histo_->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=histo_->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;

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
   //TH1F* hEoverP_corrEnergyBin[nCorrEnergyBins];

   vector<TH1F*> hEcorrOverEtrue_corrEnergyBin(nCorrEnergyBins,NULL);
   vector<TH1F*> hErawOverEtrue_corrEnergyBin(nCorrEnergyBins,NULL);
   vector<TH1F*> hPtrackOverEtrue_corrEnergyBin(nCorrEnergyBins,NULL);

   for (Int_t i = 0; i < nCorrEnergyBins; i++) {
     hEoverP_corrEnergyBin[i] = new TH1F(Form("hEoverP_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",100,0.05,2.05);
   }

   if (sampleName != "DATA") {
     for (Int_t i = 0; i < nCorrEnergyBins; i++) {
       hEcorrOverEtrue_corrEnergyBin[i] = new TH1F(Form("hEcorrOverEtrue_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",200,0.55,1.55);
       hErawOverEtrue_corrEnergyBin[i] = new TH1F(Form("hErawOverEtrue_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",200,0.55,1.55);
       hPtrackOverEtrue_corrEnergyBin[i] = new TH1F(Form("hPtrackOverEtrue_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",100,0.05,2.05);
     }
   }

   // Float_t lowMeanCut = 0.85;
   // Float_t upMeanCut = 1.15;

   // used only for MC study
   /////////////////////////
   Int_t MCtruthMatchFound = 0;
   Float_t Etrue = 0.0; 
   /////////////////////////

   // TH1F* hMeanEoverP = new TH1F("hMeanEoverP","",nCorrEnergyBins,corrEnergybinEdges.data());
   // TH1F* hModeEoverP = new TH1F("hModeEoverP","",nCorrEnergyBins,corrEnergybinEdges.data());

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


      // look for the bin in the LepGood_correctedEcalEnergy variable
      Double_t theta = 2. * atan(exp(-LepGood_eta[0])); 
      Double_t energyToUse = LepGood_correctedEcalEnergy[0]*sin(theta); 
      Int_t bin = getBinNumber(energyToUse,corrEnergybinEdges);  // this function returns negative value if bin not found

      if (bin >= 0) {

	hEoverP_corrEnergyBin[bin]->Fill(LepGood_eSuperClusterOverP[0]);

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

   }

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

void plotDistribution(const string &sampleName, const vector<Float_t> &corrEnergybinEdges, TH1F* hPeak, TH1F* hSigma, const string &hNameID, const string &dirName) {

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
    
  UInt_t nBins = corrEnergybinEdges.size() -1;

  for (UInt_t i = 0; i < nBins; i++) {

    htmp = (TH1F*)f->Get(Form("h%s_corrEnergyBin%1.0fTo%1.0f",hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1]));  
    if (!htmp) {
      cout << "Error: histogram not found in file ' " << fileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hist = (TH1F*)htmp->Clone();

    hist->SetStats(0);  
    hist->Draw("HE");
    if (sampleName == "DATA") {
      if (corrEnergybinEdges[i] > 349.9) hist->Rebin(3); 
      else if (corrEnergybinEdges[i] > 199.9) hist->Rebin(2); // when using E instead of Et 249.9 is ok
    } else if (hNameID != "EoverP") {
      if (hNameID == "PtrackOverEtrue") {
	hist->Rebin(2);
      } else {
	hist->GetXaxis()->SetRangeUser(0.85,1.15);
      }
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
    hist->Fit("gaus","WL I Q 0","",gaussEdgeL,gaussEdgeR);  // L: loglikelihood method, 0: do not plot this fit, Q: quiet mode (minimum printing)
    Double_t gaussNorm = hist->GetFunction("gaus")->GetParameter(0);
    Double_t gaussMean = hist->GetFunction("gaus")->GetParameter(1);
    //Double_t gaussMeanError = hist->GetFunction("gaus")->GetParError(1);
    Double_t gaussSigma = hist->GetFunction("gaus")->GetParameter(2);
    // now use crystal ball with right tail
    Double_t funcRangeLeft = gaussEdgeL;
    Double_t funcRangeRight = 2.0;
    if (hNameID != "EoverP") {
      if (hNameID != "PtrackOverEtrue") funcRangeLeft = 0.8;
      else funcRangeLeft = 0.2;
      funcRangeRight = gaussEdgeR;
    }
    TF1 *cb1;
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
    TFitResultPtr frp1 = hist->Fit(cb1,"WL I S Q B R","HE",funcRangeLeft,funcRangeRight);
    if (hPeak != NULL) {
      hPeak->SetBinContent(i+1, frp1->Parameter(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
      hPeak->SetBinError(i+1, frp1->ParError(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
    }
    if (hSigma != NULL) {
      hSigma->SetBinContent(i+1, frp1->Parameter(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
      hSigma->SetBinError(i+1, frp1->ParError(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
    }

    hist->SetTitle(Form("%1.0f < E_{T}[GeV] < %1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));
    
    c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.pdf",dirName.c_str(),hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));
    c->SaveAs(Form("%s%sdistribution_ET%1.0fTo%1.0f_%s.png",dirName.c_str(),hNameID.c_str(),corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));

  }

  delete c;
  delete htmp;
  delete hist;

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
  //if (0) ratioplot->GetYaxis()->SetRangeUser(0.4,1.4);
  ratioplot->SetMarkerStyle(8);  //medium dot  
  ratioplot->DrawCopy("E");
  cfit->SaveAs( (dirName + canvasName + ".pdf").c_str() );
  cfit->SaveAs( (dirName + canvasName + ".png").c_str() );


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


}

//=================================================================================

void plotFromFit(const string &dataSampleName, const string &MCSampleName, const vector<Float_t> &corrEnergybinEdges, const string & dirName) {
  
  Int_t nCorrEnergyBins = corrEnergybinEdges.size() -1;

  TH1F* hPeakEoverPdata = new TH1F("hPeakEoverPdata","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEoverPdata = new TH1F("hSigmaEoverPdata","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hPeakEoverPmc = new TH1F("hPeakEoverPmc","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEoverPmc = new TH1F("hSigmaEoverPmc","",nCorrEnergyBins,corrEnergybinEdges.data());

  plotDistribution(dataSampleName, corrEnergybinEdges, hPeakEoverPdata, hSigmaEoverPdata, "EoverP",dirName);
  plotDistribution(MCSampleName, corrEnergybinEdges, hPeakEoverPmc, hSigmaEoverPmc, "EoverP",dirName);

  drawPlotDataMC(hPeakEoverPdata, hPeakEoverPmc, MCSampleName, "corrected E_{T} [GeV]", "peak(E/P)", "modeEoverPfromFit",dirName);
  drawPlotDataMC(hSigmaEoverPdata, hSigmaEoverPmc, MCSampleName, "corrected E_{T} [GeV]", "#sigma(E/P)", "sigmaEoverPfromFit",dirName);

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

  drawPlotOnlyMC(hPeakVectorMC, legEntryName, MCSampleName, "corrected E_{T} [GeV]", "peak position", "modeMCstudy",dirName); 
  drawPlotOnlyMC(hSigmaVectorMC, legEntryName, MCSampleName, "corrected E_{T} [GeV]", "#sigma of distribution", "sigmaMCstudy",dirName); 


}

//======================================================================


Int_t main(Int_t argc, char* argv[]) {

  Int_t doAll_flag = 1;
  Int_t doLoop_flag = 1;
  Int_t doPlot_flag = 1;

  string dirName = "";   // will be a path like "plot/<name>/" . Note the ending "/" 

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
	cout << "Saving output in '" << dirName << "'" << endl;
	i++;
      }

    }

  }

  vector<string> sampleName;
  sampleName.push_back("DATA");
  sampleName.push_back("WJetsToLNu");

  vector<Float_t> corrEnergybinEdges;
  corrEnergybinEdges.push_back(25.0);
  corrEnergybinEdges.push_back(50.0);
  corrEnergybinEdges.push_back(75.0);
  corrEnergybinEdges.push_back(100.0);
  corrEnergybinEdges.push_back(150.0);
  corrEnergybinEdges.push_back(200.0);
  corrEnergybinEdges.push_back(250.0);
  //corrEnergybinEdges.push_back(300.0);
  corrEnergybinEdges.push_back(350.0);
  // corrEnergybinEdges.push_back(450.0);
  // corrEnergybinEdges.push_back(600.0);
  // corrEnergybinEdges.push_back(750.0);
  corrEnergybinEdges.push_back(900.0);

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
