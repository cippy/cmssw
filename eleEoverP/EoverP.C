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

//#define ISDATA_FLAG 1
#define DIR_TO_SAVE_PLOT "./plot/"

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

void buildChainWithFriend(TChain* chain, TChain* chFriend, string sampleName) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "DATA") {

    subSampleNameVector.push_back("SingleElectron_PromptReco_v1_runs_272021_273149");
    subSampleNameVector.push_back("SingleElectron_PromptReco_v2_runs_273150_274443");

  } else if (sampleName == "WJetsToLNu") {
    
    subSampleNameVector.push_back("WJetsToLNu_HT100to200");
    subSampleNameVector.push_back("WJetsToLNu_HT200to400");
    subSampleNameVector.push_back("WJetsToLNu_HT400to600");
    subSampleNameVector.push_back("WJetsToLNu_HT600to800");
    subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
    subSampleNameVector.push_back("WJetsToLNu_HT2500toInf");  // note, 1200to2500 bin missing                                    

  } else {

    cout << "Error: unknown sampleName " << sampleName <<". End of programme" << endl;
    exit(EXIT_FAILURE);

  }

  string treePath = "root://eoscms//eos/cms/store/cmst3/group/susy/emanuele/monox/trees/TREES_1LEPSKIM_80X/";
  
  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    string friend_treeRootFile = treePath + "friends_evVarFriend_" + subSampleNameVector[i]+ ".root";

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

void EoverP::Loop(const string sampleName, const vector<Float_t> &corrEnergybinEdges)
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
   fChain->SetBranchStatus("LepGood_eSuperClusterOverP",1);
   fChain->SetBranchStatus("LepGood_r9",1);
   fChain->SetBranchStatus("met_pt",1);

   // branches for MC study
   if (sampleName != data){
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

   string rootfileName = "EoverP_" + sampleName + ".root";

   TFile *rootFile = new TFile((rootfileName).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout << "Error: file \"" << rootfileName << "\" was not opened." << endl;
     exit(EXIT_FAILURE);
   }

   vector<TH1F*> hEoverP_corrEnergyBin(nCorrEnergyBins,NULL);
   //TH1F* hEoverP_corrEnergyBin[nCorrEnergyBins];

   for (Int_t i = 0; i < nCorrEnergyBins; i++) {
     hEoverP_corrEnergyBin[i] = new TH1F(Form("hEoverP_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]),"",100,0.05,2.05);
   }

   // Float_t lowMeanCut = 0.85;
   // Float_t upMeanCut = 1.15;

   Float_t Etrue = 0.0; // used only for MC study

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

      // here goes with the algorithm to match gen to reco electrons

      if (sampleName != "data") {

	Int_t matchFound = 0;
	Int_t i = 0;

	while (!matchFound && i < ngenLep) {

	  // start asking for e+/- coming from a W+/-
	    
	  if (fabs(genLep_pdgId[i]) == 11 && fabs(genLep_motherId[i]) == 24) {

	    TVector3 eleReco;
	    eleReco.SetPtEtaPhi(LepGood_pt[0],LepGood_eta[0],LepGood_phi[0]);
	    TVector3 eleGen;
	    eleGen.SetPtEtaPhi(genLep_pt[0],genLep_eta[0],genLep_phi[0]);

	    // now match is ok if deltaR < 0.51, which is ~ 3 crystal along eta
	    // in this case, compute Etrue and set matchFound flag as 1

	    if (eleGen.DeltaR(eleReco) < 0.051) {

	      Etrue = eleGen.Mag();  // neglect mass since we have electrons
	      matchFound = 1;  

	    }
	    
	  }

	  i++;

	}

      } 

      // end of loop to match gen and reco electrons

      Int_t binFound = 0;
      Int_t bin = 0;
      while (!binFound && bin < nCorrEnergyBins) {
	if (LepGood_correctedEcalEnergy[0] < corrEnergybinEdges[bin+1]) {
	  hEoverP_corrEnergyBin[bin]->Fill(LepGood_eSuperClusterOverP[0]);
	  binFound = 1;
	}
	bin++;
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

     // decided to store only the distribution, without aleady filling histograms with E/P or whatever at the peak

     // //Double_t stdDev = hEoverP_corrEnergyBin[i]->GetRMS();
     // //     Double_t stdDev = getEffectiveSigma(hEoverP_corrEnergyBin[i]);
     // //if (stdDev <= 0.000001) stdDev = hEoverP_corrEnergyBin[i]->GetRMS(); // if 0 (being double, better to ask < epsilon)  
     // hMeanEoverP->SetBinContent(i+1,hEoverP_corrEnergyBin[i]->GetMean());  // try also with real peak position                                    
     // //hMeanEoverP->SetBinError(i+1,stdDev);  // use RMS as uncertainty, can always get back to poisson uncertainty by taking square root
     // hMeanEoverP->SetBinError(i+1,hEoverP_corrEnergyBin[i]->GetMeanError());  // use uncertainty on the mean 
     // // again with mean in restricted range to get the peak
     // TH1F* hTmp = new TH1F(*(hEoverP_corrEnergyBin[i]));
     // hTmp->GetXaxis()->SetRangeUser(lowMeanCut,upMeanCut);     
     // hModeEoverP->SetBinContent(i+1,hTmp->GetMean());  // try also with real peak position                              
     // //hModeEoverP->SetBinError(i+1,stdDev);  // same uncertainty as above (RMS of whole distribution)
     // hModeEoverP->SetBinError(i+1,hTmp->GetMeanError());  // use uncertainty on the mean in the restricted range (so less entries and larger uncertainty) 
     // delete hTmp;

   }

   // hMeanEoverP->SetStats(0);
   // hMeanEoverP->GetXaxis()->SetTitle("corrected E [GeV]");
   // hMeanEoverP->GetXaxis()->SetTitleSize(0.06);
   // hMeanEoverP->GetXaxis()->SetTitleOffset(0.75);
   // hMeanEoverP->GetYaxis()->SetTitle("< E/P >");
   // hMeanEoverP->GetYaxis()->SetTitleSize(0.055);
   // hMeanEoverP->GetYaxis()->SetTitleOffset(0.75);

   // hModeEoverP->SetStats(0);
   // hModeEoverP->GetXaxis()->SetTitle("corrected E [GeV]");
   // hModeEoverP->GetXaxis()->SetTitleSize(0.06);
   // hModeEoverP->GetXaxis()->SetTitleOffset(0.75);
   // hModeEoverP->GetYaxis()->SetTitle("Mode of E/P ");
   // hModeEoverP->GetYaxis()->SetTitleSize(0.055);
   // hModeEoverP->GetYaxis()->SetTitleOffset(0.75);

   rootFile->Write();
   rootFile->Close();
   delete rootFile;


}

void getTexMCSampleName(const string &MCSampleName, string &texMCSampleName) {
  
  if (MCSampleName == "WJetsToLNu") texMCSampleName = "W(l#nu)+jets";

}

void plotEoverPdistribution(const string &sampleName, const vector<Float_t> &corrEnergybinEdges, TH1F* hPeakEoverP, TH1F* hSigmaEoverP) {

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 
  TVirtualFitter::SetDefaultFitter("Minuit");

  string fileName = "EoverP_" + sampleName + ".root";

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

    htmp = (TH1F*)f->Get(Form("hEoverP_corrEnergyBin%1.0fTo%1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));  
    if (!htmp) {
      cout << "Error: histogram not found in file ' " << fileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hist = (TH1F*)htmp->Clone();

    hist->SetStats(0);  
    hist->Draw("HE");
    if (sampleName == "DATA") {
      if (corrEnergybinEdges[i] > 349.9) hist->Rebin(3); 
      else if (corrEnergybinEdges[i] > 249.9) hist->Rebin(2);
    }

    // fitting:
    // do a first fit with a simple gaussian in the core
    Double_t gaussEdgeL = 0.9;  //left side of the gaussian to be used in the fit (I use a variable so that I change this value only once)
    Double_t gaussEdgeR = 1.1;  //right side ...
    if (corrEnergybinEdges[i] > 349.9) {
      gaussEdgeL = 0.8;
      gaussEdgeR = 1.2;
    }
    hist->Fit("gaus","L 0 Q","",gaussEdgeL,gaussEdgeR);  // L: loglikelihood method, 0: do not plot this fit, Q: quiet mode (minimum printing)
    Double_t gaussNorm = hist->GetFunction("gaus")->GetParameter(0);
    Double_t gaussMean = hist->GetFunction("gaus")->GetParameter(1);
    //Double_t gaussMeanError = hist->GetFunction("gaus")->GetParError(1);
    Double_t gaussSigma = hist->GetFunction("gaus")->GetParameter(2);
    // now use crystal ball with right tail
    Double_t funcRangeLeft = gaussEdgeL;
    Double_t funcRangeRight = 2.0;
    TF1 *cb1 = new TF1("cb1",&myCrystalBallRightTail,funcRangeLeft,funcRangeRight,5);  // last parameter is the number of free parameters
    cb1->SetParNames("alpha","n","mu","sigma","N");  
    cb1->SetParLimits(cb1->GetParNumber("n"),0.1,15); 
    cb1->SetParLimits(cb1->GetParNumber("alpha"),0.01,10);
    cb1->SetParameters((gaussEdgeR-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm);
    TFitResultPtr frp1 = hist->Fit(cb1,"WL I S Q B R","HE",funcRangeLeft,funcRangeRight);
    hPeakEoverP->SetBinContent(i+1, frp1->Parameter(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
    hPeakEoverP->SetBinError(i+1, frp1->ParError(2)); // 2 is mu (starts with alpha, which is parameter number 0)  
    hSigmaEoverP->SetBinContent(i+1, frp1->Parameter(3)); // 2 is mu (starts with alpha, which is parameter number 0)  
    hSigmaEoverP->SetBinError(i+1, frp1->ParError(3)); // 2 is mu (starts with alpha, which is parameter number 0)  

    hist->SetTitle(Form("%1.0f < E[GeV] < %1.0f",corrEnergybinEdges[i],corrEnergybinEdges[i+1]));

    c->SaveAs(Form("%sEoverPdistribution_E%1.0fTo%1.0f_%s.pdf",DIR_TO_SAVE_PLOT,corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));
    c->SaveAs(Form("%sEoverPdistribution_E%1.0fTo%1.0f_%s.png",DIR_TO_SAVE_PLOT,corrEnergybinEdges[i],corrEnergybinEdges[i+1],sampleName.c_str()));

  }

  delete c;
  delete htmp;
  delete hist;

}

//==========================================================================================


void drawPlot(TH1F* hdata, TH1F* hmc, const string& MCSampleName, const string& xAxisName, const string& yAxisName, const string& canvasName) {

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
  ratioplot->SetMarkerStyle(8);  //medium dot  
  ratioplot->DrawCopy("E");
  string dir = string(DIR_TO_SAVE_PLOT);
  cfit->SaveAs( (dir + canvasName + ".pdf").c_str() );
  cfit->SaveAs( (dir + canvasName + ".png").c_str() );


}


//=================================================================================

void plotEoverPfromFit(const string &dataSampleName, const string &MCSampleName, const vector<Float_t> &corrEnergybinEdges) {
  
  Int_t nCorrEnergyBins = corrEnergybinEdges.size() -1;

  TH1F* hPeakEoverPdata = new TH1F("hPeakEoverPdata","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEoverPdata = new TH1F("hSigmaEoverPdata","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hPeakEoverPmc = new TH1F("hPeakEoverPmc","",nCorrEnergyBins,corrEnergybinEdges.data());
  TH1F* hSigmaEoverPmc = new TH1F("hSigmaEoverPmc","",nCorrEnergyBins,corrEnergybinEdges.data());

  plotEoverPdistribution(dataSampleName, corrEnergybinEdges, hPeakEoverPdata, hSigmaEoverPdata);
  plotEoverPdistribution(MCSampleName, corrEnergybinEdges, hPeakEoverPmc, hSigmaEoverPmc);

  drawPlot(hPeakEoverPdata, hPeakEoverPmc, MCSampleName, "corrected E [GeV]", "mode of E/P", "modeEoverPfromFit");
  drawPlot(hSigmaEoverPdata, hSigmaEoverPmc, MCSampleName, "corrected E [GeV]", "#sigma(E/P)", "sigmaEoverPfromFit");

}

//======================================================================


void plotEoverP(const string &dataSampleName, const string &MCSampleName, const string &histName, const string &yAxisName) {

  // currently not using this function any more


  // data(MC)FileName --> name for data (MC) root file
  // histName --> name of histogram inside root file
  // yAxisName --> name to assign to y axis (somehow related to histName)


  string dataFileName = "EoverP_" + dataSampleName + ".root";
  string MCFileName = "EoverP_" + MCSampleName + ".root";

  TCanvas *c = new TCanvas("c","",700,700);  

  TH1F* htmp = NULL;
  TH1F* hdata = NULL;
  TH1F* hmc = NULL;
    
  TFile* f = TFile::Open(dataFileName.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<dataFileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  htmp = (TH1F*)f->Get(histName.c_str());
    
  if (!htmp) {
    cout << "Error: histogram not found in file ' " << dataFileName << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  hdata = (TH1F*)htmp->Clone();

  TFile* fmc = TFile::Open(MCFileName.c_str(),"READ");
  if (!fmc || !fmc->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<MCFileName<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
    
  htmp = (TH1F*)fmc->Get(histName.c_str());
    
  if (!htmp) {
    cout << "Error: histogram not found in file ' " << MCFileName << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  hmc = (TH1F*)htmp->Clone();
    

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
  ratioplot->GetXaxis()->SetTitle("corrected E [GeV]");
  ratioplot->GetXaxis()->SetTitleSize(0.14);
  ratioplot->GetXaxis()->SetTitleOffset(0.8);
  ratioplot->GetYaxis()->SetLabelSize(0.10);
  ratioplot->GetYaxis()->SetTitle("data / MC");
  ratioplot->GetYaxis()->SetTitleSize(0.12);
  ratioplot->GetYaxis()->SetTitleOffset(0.45);
  ratioplot->GetYaxis()->CenterTitle();
  ratioplot->GetYaxis()->SetNdivisions(011);
  ratioplot->SetMarkerStyle(8);  //medium dot  
  ratioplot->DrawCopy("E");
  string dir = string(DIR_TO_SAVE_PLOT);
  c->SaveAs( (dir + histName + ".pdf").c_str() );
  c->SaveAs( (dir + histName + ".png").c_str() );

}


Int_t main(Int_t argc, char* argv[]) {

  Int_t doAll_flag = 1;
  Int_t doLoop_flag = 1;

  if (argc > 1) {

    for (Int_t i = 1; i < argc; i++) {

      string thisArgument(argv[i]);
      if (thisArgument == "-nl") {
	cout << "Passing option -nl: skip Loop on ntuples" << endl;
	doLoop_flag = 0;  // -nl --> no Loop
	doAll_flag = 0;
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
      tree.Loop(sampleName[i], corrEnergybinEdges);

      delete chain;
      delete chFriend;

    }

  }

  plotEoverPfromFit(sampleName[0],sampleName[1], corrEnergybinEdges);

  // plotEoverP(sampleName[0],sampleName[1],"hMeanEoverP","< E/P >");
  // plotEoverP(sampleName[0],sampleName[1],"hModeEoverP"," mode of E/P");
  // plotEoverPdistribution(sampleName[0], corrEnergybinEdges);
  // plotEoverPdistribution(sampleName[1], corrEnergybinEdges);

  return 0;

}
