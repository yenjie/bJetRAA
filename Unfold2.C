#include <iostream>
#include <stdio.h>

#include <TRandom2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "utilities.h"
#include "bayesianUnfold.h"
#include "prior.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
//==============================================================================

void Unfold2(int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 0,char *spectraFileName = (char*)"pbpb_spectra_akPu3PF.root",double recoJetPtCut = 60.,double trackMaxPtCut = 0, int nBayesianIter = 4, int doBjets=1, int doTrigEffCorr=1) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
{
  
  gStyle->SetErrorX(0.);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);	
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  // input files
  char *fileNamePP_mc = NULL;
  char *fileNamePbPb_mc = NULL;
  char *fileNamePP_data = NULL;
  char *fileNamePbPb_data = NULL;
  

  // pp file needs replacing
  if(doBjets)fileNamePP_data = (char*)"~/Work/bTagging/outputTowardsFinal/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root";
  //else fileNamePP_data = (char*)"../spectra/rawIncSpectra_MB_fixedBinSize_v4.root";
  else fileNamePP_data = (char*)"../spectra/rawIncSpectra_MB_varBinSize_v4.root";
  //if(doBjets)fileNamePbPb_data = (char*)"~/Work/bTagging/outputTowardsFinal/AltBinningV6_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root";
  if(doBjets)fileNamePbPb_data = (char*)"~/Work/bTagging/outputTowardsFinal/bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root";
  //else fileNamePbPb_data = (char*)"../spectra/rawIncSpectra_MB_fixedBinSize_v4.root";
  else fileNamePbPb_data = (char*)"../spectra/rawIncSpectra_MB_varBinSize_v4.root";
  if(doBjets) fileNamePP_mc = (char*)"~/Work/bTagging/histos/ppMC_ppReco_ak3PF_BjetTrig_noIPupperCut.root";
  else fileNamePP_mc = (char*)"~/Work/bTagging/histos/ppMC_ppReco_ak3PF_QCDjetTrig_noIPupperCut.root";
  if(doBjets)fileNamePbPb_mc = (char*) "~/Work/bTagging/histos/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root";
  else fileNamePbPb_mc = (char*) "~/Work/bTagging/histos/PbPbQCDMC_pt30by3_ipHICalibCentWeight_noTrig.root";

  // grab ntuples
  TFile *infPbPb_mc = new TFile(fileNamePbPb_mc);
  TFile *infPP_mc = new TFile(fileNamePP_mc);


  string bJetString = "Inc";
  if(doBjets) bJetString = "bJets";

  // Output file
  TFile *pbpb_Unfo;
  if (isMC) pbpb_Unfo = new TFile(Form("pbpb_Unfo_%s_MC_v2.root",algoName[algo]),"RECREATE");
  else pbpb_Unfo  = new TFile(Form("pbpb_Unfo_%s_jtpt%.0f_%s_v4.root",algoName[algo],recoJetPtCut,bJetString.c_str()),"RECREATE");

  // Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];
		
  // Initialize Histograms   
	
  for (int i=0;i<=nbins_cent;i++) uhist[i] = new UnfoldingHistos(i);
	
  // Initialize reweighting functions
  /*
  TCut dataSelection;
  TCut dataSelectionPP;
  TCut TriggerSelectionPP;
  TCut TriggerSelectionPbPb80;

  if(doBjets)dataSelection = "weight*(abs(refparton_flavorForB)==5&&abs(jteta)<2)";
  else dataSelection = "weight*(abs(jteta)<2)";
  */

  if (isMC) cout<<"This is a MC closure test"<<endl;
  else cout<< "This is a data analysis"<<endl;    		     
	     	
  // Setup jet data branches, basically the jet tree branches are assigned to this object when we loop over the events
	
  JetDataPbPb *dataPbPb   = new JetDataPbPb(fileNamePbPb_mc,(char*)"nt"); // PbPb data	
  JetDataPP *dataPP = new JetDataPP(fileNamePP_mc,(char*)"nt");	// pp data
	
  TFile *fSpectra(0);		
  if (useSpectraFromFile||useMatrixFromFile){
    fSpectra = new TFile(spectraFileName,"read");
  }
  
  // Come back to the output file dir
  pbpb_Unfo->cd();

  // Get Jet spectra from data file
  cout <<"Reading data..."<<endl;

  // This doesn't seem to be relevant for the moment -Matt
  /*	

  TTree *tPbPbJet = (TTree*)infPbPb_mc->Get("nt");
  TTree *tPPJet  = (TTree*)infPP_mc->Get("nt");


  TCanvas * cInput = new TCanvas("cInput","Input",800,400);
  cInput->Divide(2,1);
		
  cout <<"Spectra..."<<endl;	
	
  for (int i=0;i<=nbins_cent;i++){
    cout <<nbins_cent<<endl;
    TCut centCut = Form("bin<%.0f&&bin>=%.0f",boundaries_cent[i+1],boundaries_cent[i]);
    if (useSpectraFromFile) {
      uhist[i]->hMeas = (TH1F*)fSpectra->Get(Form("hMeas_cent%d",i));
    } else {
      if (!isMC) {
	tPbPbJet->Project(Form("hMeas_cent%d",i),"jtptB", dataSelection&&centCut&&TriggerSelectionPbPb80);
      }   
    }
		
    if (useMatrixFromFile) {
      cout <<"Matrix"<<endl;
      uhist[i]->hMatrixFit = (TH2F*) fSpectra->Get(Form("hMatrixFit_cent%d",i));
      uhist[i]->hMeasMatch = (TH1F*)((TH2F*) fSpectra->Get(Form("hMatrixFit_cent%d",i)))->ProjectionY();
      uhist[i]->hMeasMatch->Divide(uhist[i]->hMeas);
    } else {
      uhist[i]->hMeasMatch = 0;
    }
    uhist[i]->hMeas->Draw();
  }

  if (!isMC) tPPJet->Project(Form("hMeas_cent%d",nbins_cent),"jtpt",dataSelectionPP&&TriggerSelectionPP);
  */
  cout <<"MC..."<<endl;	
  
  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
 

  // if you change the binning you have to change these, too
  // inclusive trig eff
  /*
    float trigEffInc[6]={0.777265,
			 0.95765,
			 0.998357,
			 0.999941,
			 1.,
			 1.};
  */

    // b-jet trig eff
    /*
    float trigEffbJet[6]={0.660276,
		      0.908997,
		      0.980793,
		      0.998767,
		      0.999442,
		      1.};
    */



 
		
  // Fill PbPb MC   
  if (!useMatrixFromFile) {
    for (Long64_t jentry2=0; jentry2<dataPbPb->tJet->GetEntries();jentry2++) {
      dataPbPb->tJet->GetEntry(jentry2);
      
      // change when we switch to centrality binning
      int cBin = 0;
      
      //int cBin = hCent->FindBin(dataPbPb->bin)-1;
      /*
	if (cBin>=nbins_cent) continue;
	if (cBin==-1) continue;
      */
      
      if ( dataPbPb->refpt  < 0. ) continue;
      if ( dataPbPb->jteta  > 2. || dataPbPb->jteta < -2. ) continue;
      if ( dataPbPb->refpt<0) dataPbPb->refpt=0;
      if (doBjets && fabs(dataPbPb->refparton_flavorForB)!=5) continue;
      //if (doBjets&& dataPbPb->discr_ssvHighEff<2) continue;
      if (doBjets && dataPbPb->jtptB < recoJetPtCut) continue;
      if (!doBjets && dataPbPb->jtptA < recoJetPtCut) continue;
      //if (!doTrigEffCorr && dataPbPb->isTrig <1) continue;
      if ( dataPbPb->isTrig <1) continue;
      
      //if(!doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptA>120) continue;
      //if(doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptB>120) continue;

      float weight = dataPbPb->weight;

      if(doTrigEffCorr){
	for(int iBin = 0; iBin<nbins_rec; iBin++){
	  float myJtPt = 0.;
	  if(doBjets) myJtPt = dataPbPb->jtptB;
	  else myJtPt = dataPbPb->jtptA;
	  if(myJtPt > boundaries_rec[iBin] && myJtPt < boundaries_rec[iBin+1]){
	    if(doBjets) weight/= trigEffbJet[iBin];
	    else weight/= trigEffInc[iBin];
	  }							  
	}
      }
      if (!isMC||jentry2 % 2 == 1) {
	if(doBjets)uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtptB,weight);
	else uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtptA,weight);
      }	  
      if (jentry2 % 2 == 0) {
	uhist[cBin]-> hGen->Fill(dataPbPb->refpt,weight);   
	if(doBjets)uhist[cBin]-> hMeas->Fill(dataPbPb->jtptB,weight);  	 
	else uhist[cBin]-> hMeas->Fill(dataPbPb->jtptA,weight);  	 
	//uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-i)),weight); 
	// FIXME!!!!!!  i is supposed to be a loop over centrality !!!!
	if(doBjets)uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtptB*(1.+0.02/nbins_cent*(nbins_cent-0)),weight); 
	else uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtptA*(1.+0.02/nbins_cent*(nbins_cent-0)),weight); 
      }
    }

    //pp will just fill the last index of the centrality array

    // fill pp MC
    for (Long64_t jentry2=0; jentry2<dataPP->tJet->GetEntries();jentry2++) {
      dataPP->tJet->GetEntry(jentry2);
      
      if ( dataPP->refpt<0) continue;
      if ( dataPP->jteta  > 2. || dataPP->jteta < -2. ) continue;
      if ( dataPP->refpt<0) dataPP->refpt=0;
      if ( doBjets && fabs(dataPP->refparton_flavorForB)!=5) continue;
      //if ( doBjets && dataPP->discr_ssvHighEff<2) continue;
      if ( dataPP->jtpt < recoJetPtCut) continue;
      
      if (!isMC||jentry2 % 2 == 1) {
	uhist[nbins_cent]-> hMatrix->Fill(dataPP->refpt,dataPP->jtpt,dataPP->weight);
      }	  
      if (jentry2 % 2 == 0) {
	uhist[nbins_cent]-> hGen->Fill(dataPP->refpt,dataPP->weight);   
	uhist[nbins_cent]-> hMeas->Fill(dataPP->jtpt,dataPP->weight); 
      }           
    }
  }
  /*
  for (int i=0;i<=nbins_cent;i++){
    for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	float binContent = uhist[i]->hGen->GetBinContent(x);
	float binError = uhist[i]->hGen->GetBinError(x);
	float binWidth =  uhist[i]->hGen->GetXaxis()->GetBinWidth(x);
	uhist[i]->hGen->SetBinContent(x,binContent/binWidth);
	uhist[i]->hGen->SetBinError(x,binError/binWidth);
    }
  }

  for (int i=0;i<=nbins_cent;i++){
    for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	float binContent = uhist[i]->hMeas->GetBinContent(x);
	float binError = uhist[i]->hMeas->GetBinError(x);
	float binWidth =  uhist[i]->hMeas->GetXaxis()->GetBinWidth(x);
	uhist[i]->hMeas->SetBinContent(x,binContent/binWidth);
	uhist[i]->hMeas->SetBinError(x,binError/binWidth);
    }
  }

  for (int i=0;i<=nbins_cent;i++){
    for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
      for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
	float binContent = uhist[i]->hMatrix->GetBinContent(x,y);
	float binError = uhist[i]->hMatrix->GetBinError(x,y);
	float binWidth2 =  uhist[i]->hMatrix->GetXaxis()->GetBinWidth(x)*uhist[i]->hMatrix->GetYaxis()->GetBinWidth(y);
	uhist[i]->hMatrix->SetBinContent(x,y,binContent/binWidth2);
	uhist[i]->hMatrix->SetBinError(x,y,binError/binWidth2);
      }      
    }
  }	
  */


  if (isMC==0) {
    // Use measured histogram from Matt & Kurt's file
	   
    // PbPb file:

    TFile *infMatt = new TFile(fileNamePbPb_data);
    TH1F *hMattPbPb = NULL;
    TH1F *hTagEffPbPb = NULL;

    if(doBjets){
      hMattPbPb = (TH1F*) infMatt->Get("hRawBData");
      hTagEffPbPb = (TH1F*) infMatt->Get("hBEfficiencyMC");
    }
    else hMattPbPb = (TH1F*) infMatt->Get("hPbPb");
    //divideBinWidth(hMattPbPb);


    
    for (int i=1;i<=hMattPbPb->GetNbinsX();i++)
      {
     	float binContent = hMattPbPb->GetBinContent(i);  
	float binError = hMattPbPb->GetBinError(i); 

	if(doBjets){
	  float tagEff =hTagEffPbPb->GetBinContent(i);
	  float tagEffErr =     hTagEffPbPb->GetBinError(i);   
	  
	  if(tagEff>0){
	    // careful of the order here!
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	    binContent /= tagEff;
	  }
	  else cout<<" TAGEFF = 0"<<endl;	  
	}

	float binCenter = hMattPbPb->GetBinCenter(i);  
	if(binCenter - hMattPbPb->GetBinWidth(i)/2.  < recoJetPtCut) continue;
	
	int ibin=0;
	
	for(ibin=1;ibin<=uhist[0]->hMeas->GetNbinsX();ibin++){
	  float testLowEdge = uhist[0]->hMeas->GetBinLowEdge(ibin);
	  float testBinWidth = uhist[0]->hMeas->GetBinWidth(ibin);
	  if(binCenter>testLowEdge && binCenter < testLowEdge+testBinWidth) break;
	}
	
	
	if(doTrigEffCorr){
	  float trigEff = 0;
	  if(doBjets) trigEff = trigEffbJet[ibin-1];
	  else  trigEff = trigEffInc[ibin-1];

	  cout<<" binCenter = "<<binCenter<<" trigEff "<<trigEff<<endl;

	  if(trigEff>0){
	    // careful of the order here!
	    binContent /= trigEff;
	    binError /= trigEff;
	  }
	  else cout<<" TRIGEFF = 0"<<endl;	  
	}


     
	uhist[0]->hMeas->SetBinContent(ibin,binContent);  
	uhist[0]->hMeas->SetBinError(ibin,binError);  

      }

    // pp file:
    TFile *infMattPP = new TFile(fileNamePP_data);
    TH1F *hMattPP = NULL;
    TH1F *hTagEffPP = NULL;
    if(doBjets){
      hMattPP = (TH1F*) infMattPP->Get("hRawBData");
      hTagEffPP = (TH1F*) infMattPP->Get("hBEfficiencyMC");
    }
    else hMattPP = (TH1F*) infMattPP->Get("hpp");
    //divideBinWidth(hMattPP);
	   
    for (int i=1;i<=hMattPP->GetNbinsX();i++)
      {
      	float binContent = hMattPP->GetBinContent(i);  
	float binError = hMattPP->GetBinError(i);  
	
	if(doBjets){
	  float tagEff =hTagEffPP->GetBinContent(i);
	  float tagEffErr =     hTagEffPP->GetBinError(i);   
	  if(tagEff>0){
	    // careful of the order here!
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	    binContent /= tagEff;
	  }
	  else cout<<" TAGEFF = 0"<<endl;	  
	}
	
     	float binCenter = hMattPP->GetBinCenter(i);  
	if(binCenter - hMattPP->GetBinWidth(i)/2.  < recoJetPtCut) continue;

	int ibin=0;

	for(ibin=1;ibin<=uhist[nbins_cent]->hMeas->GetNbinsX();ibin++){
	  float testLowEdge = uhist[nbins_cent]->hMeas->GetBinLowEdge(ibin);
	  float testBinWidth = uhist[nbins_cent]->hMeas->GetBinWidth(ibin);
	  if(binCenter>testLowEdge && binCenter < testLowEdge+testBinWidth) break;
	}
	uhist[nbins_cent]->hMeas->SetBinContent(ibin,binContent);  
	uhist[nbins_cent]->hMeas->SetBinError(ibin,binError);  

    
	/*
	cout<<" i "<<i<<endl;
	int newBin = i+uhist[nbins_cent]->hMeas->FindBin(61)-1;

	cout<<" newBin "<<newBin<<endl;
	cout<<"hMattPP->GetBinCenter(i) "<<hMattPP->GetBinCenter(i)<<endl;
	cout<<"uhist[nbins_cent]->hMeas->GetBinCenter(newBin) "<<uhist[nbins_cent]->hMeas->GetBinCenter(newBin)<<endl;
	*/
      }

  }



  cout <<"Response Matrix..."<<endl;
	
  TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",800,400);
  cMatrix->Divide(2,1);

  for (int i=0;i<=nbins_cent;i++){
    cMatrix->cd(i+1);
    if (!useMatrixFromFile) {
      TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
      f->SetParameters(1e10,-8.8,40);
      for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
	double sum=0;
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	  if (uhist[i]->hMatrix->GetBinContent(x,y)<=1*uhist[i]->hMatrix->GetBinError(x,y)) {
	    uhist[i]->hMatrix->SetBinContent(x,y,0);
	    uhist[i]->hMatrix->SetBinError(x,y,0);
	  }
	  sum+=uhist[i]->hMatrix->GetBinContent(x,y);
	}
				
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {	   
	  double ratio = 1;
	  uhist[i]->hMatrix->SetBinContent(x,y,uhist[i]->hMatrix->GetBinContent(x,y)*ratio);
	  uhist[i]->hMatrix->SetBinError(x,y,uhist[i]->hMatrix->GetBinError(x,y)*ratio);
	}
      }
    }
    uhist[i]->hResponse = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponse_cent%d",i));
    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {
	if (uhist[i]->hResponse->GetBinContent(x,y)<=0*uhist[i]->hResponse->GetBinError(x,y)) {
	  uhist[i]->hResponse->SetBinContent(x,y,0);
	  uhist[i]->hResponse->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponse->GetBinContent(x,y);
      }
			
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio = uhist[i]->hMeas->GetBinContent(y)/sum;
	if (uhist[i]->hMeas->GetBinContent(y)==0) ratio = 1e-100/sum;
      }
    }
		
    uhist[i]->hResponseNorm = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponseNorm_cent%d",i));
    for (int x=1;x<=uhist[i]->hResponseNorm->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {
	if (uhist[i]->hResponseNorm->GetBinContent(x,y)<=0*uhist[i]->hResponseNorm->GetBinError(x,y)) {
	  uhist[i]->hResponseNorm->SetBinContent(x,y,0);
	  uhist[i]->hResponseNorm->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponseNorm->GetBinContent(x,y);
      }
			
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	uhist[i]->hResponseNorm->SetBinContent(x,y,uhist[i]->hResponseNorm->GetBinContent(x,y)*ratio);
	uhist[i]->hResponseNorm->SetBinError(x,y,uhist[i]->hResponseNorm->GetBinError(x,y)*ratio);
      }
    }
		
    uhist[i]->hResponse->Draw("colz");
		
    if (!useMatrixFromFile) uhist[i]->hMatrixFit = uhist[i]->hMatrix;
    uhist[i]->hMatrixFit->SetName(Form("hMatrixFit_cent%d",i));
  }


 
  pbpb_Unfo->cd();
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  //char chmet[100]; 
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
  TCanvas * cPbPb = new TCanvas("cPbPb","Comparison",1200,600);
  cPbPb->Divide(2,1); 
  cPbPb->cd(1);
	
  TH1F *hRecoBW[nbins_cent+1], *hRecoBinByBinBW[nbins_cent+1], *hMeasBW[nbins_cent+1], *hGenBW[nbins_cent+1];


  for (int i=0;i<=nbins_cent;i++) {
    cPbPb->cd(i+1)->SetLogy();   
    // Do Bin-by-bin
    TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); 
    TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",90,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    delete hBinByBinCorRaw;
    delete hMCGen;
    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
		
    // Do unfolding
    //if (isMC) uhist[i]->hMeas = (TH1F*)uhist[i]->hMatrix->ProjectionY()->Clone(Form("hMeas_cent%d",i));
    prior myPrior(uhist[i]->hMatrixFit,uhist[i]->hMeas,0);
    myPrior.unfold(uhist[i]->hMeas,1);
    TH1F *hPrior;//=(TH1F*) functionHist(fPow,uhist[i]->hMeas,Form("hPrior_cent%d",i));
    hPrior = (TH1F*)uhist[i]->hGen->Clone("hPrior");//(TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));
    removeZero(hPrior);
		
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    bayesianUnfold myUnfolding(uhist[i]->hMatrixFit,myPrior.hPrior,0);
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;

    // Iteration Systematics
    for (int j=2;j<=7;j++)
      {

	bayesianUnfold myUnfoldingSys(uhist[i]->hMatrixFit,hPrior,0);
	myUnfoldingSys.unfold(uhist[i]->hMeas,j);
	uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
      }
    
    
    uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));
    
    if (doToy) {
      TCanvas *cToy = new TCanvas("cToy","toy",600,600);
      cToy->cd();
      int nExp=1000;
      TH1F *hTmp[nbins_truth+1];
      TH1F *hTmp2[nbins_truth+1];
      for (int j=1;j<=nbins_truth;j++) {
	hTmp[j] = new TH1F(Form("hTmp%d",j),"",200,0,10.+uhist[i]->hReco->GetBinContent(j)*2);
	hTmp2[j] = new TH1F(Form("hTmp2%d",j),"",200,0,10.+uhist[i]->hRecoBinByBin->GetBinContent(j)*2);
      }
      for (int exp =0; exp<nExp; exp++) {
	TH1F *hToy = (TH1F*)uhist[i]->hMeas->Clone();   
	TH2F *hMatrixToy = (TH2F*)uhist[i]->hMatrixFit->Clone();
	hToy->SetName("hToy");
	if (exp%100==0) cout <<"Pseudo-experiment "<<exp<<endl;
	for (int j=1;j<=hToy->GetNbinsX();j++) {
	  double value = gRandom->Poisson(uhist[i]->hMeas->GetBinContent(j));
	  hToy->SetBinContent(j,value);
	}
				
	for (int j=1;j<=hMatrixToy->GetNbinsX();j++) {
	  for (int k=1;k<=hMatrixToy->GetNbinsY();k++) {
	    double value = gRandom->Gaus(uhist[i]->hMatrixFit->GetBinContent(j,k),uhist[i]->hMatrixFit->GetBinError(j,k));
	    hMatrixToy->SetBinContent(j,k,value);
	  }
	}

	prior myPriorToy(hMatrixToy,hToy,0.0);
	myPriorToy.unfold(hToy,1);
	bayesianUnfold myUnfoldingToy(hMatrixToy,myPriorToy.hPrior,0.0);
	myUnfoldingToy.unfold(hToy,nBayesianIter);
	TH1F *hRecoTmp = (TH1F*) myUnfoldingToy.hPrior->Clone();
				
	for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
	  hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
	}
	delete hToy;
	delete hRecoTmp;
	delete hMatrixToy;
      }
      TF1 *fGaus = new TF1("fGaus","[0]*TMath::Gaus(x,[1],[2])");
      for (int j=1;j<=nbins_truth;j++)
	{

	  f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
	  if (hTmp[j]->GetMean()>0) {
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
	  }	       
	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
	  if (hTmp2[j]->GetMean()>0) {
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hRecoBinByBin->SetBinError(j,f->GetParameter(2));
	  }	       
	  delete hTmp[j];
	  delete hTmp2[j];
	}
      cPbPb->cd(i+1);
    }

    uhist[i]->hMeas->SetMarkerStyle(20);
    uhist[i]->hMeas->SetMarkerColor(2);
    uhist[i]->hReco->SetMarkerStyle(25);
    uhist[i]->hReco->SetName(Form("hReco_cent%d",i));
    
    uhist[i]->hReco->SetXTitle("p_{T} (GeV/c)");    
    uhist[i]->hReco->SetYTitle("Counts");    
    uhist[i]->hReco->GetXaxis()->SetNdivisions(505);
    //uhist[i]->hReco->Draw("");    
    uhist[i]->hGen->SetLineWidth(2);
    uhist[i]->hGen->SetLineColor(2);
    //if(isMC)uhist[i]->hGen->Draw("hist same");
    //uhist[i]->hReco->Draw("same");    
    uhist[i]->hRecoBinByBin->SetMarkerStyle(28);
    uhist[i]->hRecoBinByBin->Draw("same");    
    uhist[i]->hReco->SetAxisRange(recoJetPtCut,300);
    TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
    hReproduced->SetMarkerColor(4);
    hReproduced->SetMarkerStyle(24);
    //uhist[i]->hMeas->Draw("same");    

    hRecoBW[i] = (TH1F*)uhist[i]->hReco->Clone(Form("hReco%d",i));
    hRecoBinByBinBW[i] = (TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hRecoBinByBin%d",i));
    hMeasBW[i] = (TH1F*)uhist[i]->hMeas->Clone(Form("hMeas%d",i));
    if(isMC)hGenBW[i] = (TH1F*)uhist[i]->hGen->Clone(Form("hGen%d",i));

    divideBinWidth(hRecoBW[i]);    
    if(isMC)divideBinWidth(hGenBW[i]);    
    divideBinWidth(hRecoBinByBinBW[i]);    
    divideBinWidth(hMeasBW[i]);    

    hRecoBW[i]->Draw();
    if(isMC)hGenBW[i]->Draw("hist,same");
    hRecoBinByBinBW[i]->Draw("same");
    hMeasBW[i]->Draw("same");
    
    uhist[i]->hReco->SetTitle("Baysian Unfolded");
    uhist[i]->hRecoBinByBin->SetTitle("Bin-by-bin Unfolded");

    TLegend *leg = new TLegend(0.45,0.65,0.85,0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(uhist[i]->hMeas,"Measured","pl");
    leg->AddEntry(uhist[i]->hReco,"Bayesian unfolded","pl");
    leg->AddEntry(uhist[i]->hRecoBinByBin,"Bin-by-bin unfolded","pl");
    if(isMC)leg->AddEntry(uhist[i]->hGen,"Generator level truth","l");
    leg->Draw();
  }	     
  
 
  pbpb_Unfo->Write();

  SysData systematics;
  TLine *line = new TLine(60,1,250,1);

  // Iteration systematics
  TCanvas *cIterSys = new TCanvas("cIterSys","cIterSys",1200,600);
  cIterSys->Divide(2,1);
  cIterSys->cd(2);
  TH1F *hRecoIterSysPP[100];
  TH1F *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, (char*)"hRebinPP_tmp");
  TLegend *legBayesianIterPP = myLegend(0.4,0.7,0.9,0.9);
  legBayesianIterPP->AddEntry("","PP","");
         
  for (int j=2;j<7;j++) {
    hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
    hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPP[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
      hRecoIterSysPP[j]->SetTitleOffset(1.4,"Y");
      hRecoIterSysPP[j]->SetTitleOffset(1.2,"X");
      hRecoIterSysPP[j]->SetAxisRange(0.5,1.5,"Y");
      hRecoIterSysPP[j]->Draw(); 
    } else {
      hRecoIterSysPP[j]->Draw("same");
    }
         
    checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoIterSysPP[j],0,1.1);
    legBayesianIterPP->AddEntry(hRecoIterSysPP[j],Form("Iteration %d",j),"pl");     
  }

  legBayesianIterPP->Draw();
  line->Draw();
  drawEnvelope(systematics.hSysIter[nbins_cent],(char*)"hist same");


  cIterSys->cd(1);
  TH1F *hRecoIterSysPbPb[100];
  TH1F *hRebinPbPb_tmp         = rebin(uhist[0]->hReco, (char*)"hRebinPbPb_tmp");
  TLegend *legBayesianIterPbPb = myLegend(0.4,0.7,0.9,0.9);
  legBayesianIterPbPb->AddEntry("","PbPb","");
  for (int j=2;j<7;j++) {
    hRecoIterSysPbPb[j] = rebin(uhist[0]->hRecoIterSys[j],Form("hRecoIterSysPbPb_IterSys%d",j));
    hRecoIterSysPbPb[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPbPb[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPbPb[j]->Divide(hRebinPbPb_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPbPb[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
      hRecoIterSysPbPb[j]->SetTitleOffset(1.4,"Y");
      hRecoIterSysPbPb[j]->SetTitleOffset(1.2,"X");
      hRecoIterSysPbPb[j]->SetAxisRange(0.5,1.5,"Y");
      hRecoIterSysPbPb[j]->Draw(); 
    } else {
      hRecoIterSysPbPb[j]->Draw("same");
    }
         
    checkMaximumSys(systematics.hSysIter[0],hRecoIterSysPbPb[j],0,1.1);
    legBayesianIterPbPb->AddEntry(hRecoIterSysPbPb[j],Form("Iteration %d",j),"pl");     
  }
  legBayesianIterPbPb->Draw();
  line->Draw();
  drawEnvelope(systematics.hSysIter[0],(char*)"hist same");
}





