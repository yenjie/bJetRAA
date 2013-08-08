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

void Unfold2(int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 0,char *spectraFileName = (char*)"pbpb_spectra_akPu3PF.root",double recoJetPtCut = 60,double trackMaxPtCut = 0, int nBayesianIter = 4, int doBjets=0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
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

  if(doBjets)fileNamePP_data = (char*)"/d102/mnguyen/bTaggingOutput/histos/NewFormatV4_bFractionMCTemplate_pppp1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root";
  else fileNamePP_data = (char*)"/d102/mnguyen/bTaggingOutput/histos/rawIncSpectra_MB.root";
  if(doBjets)fileNamePbPb_data = (char*) "/d102/mnguyen/bTaggingOutput/histos/AltBinningV6_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root";
  else fileNamePbPb_data = (char*)"/d102/mnguyen/bTaggingOutput/histos/rawIncSpectra_MB.root";
  if(doBjets) fileNamePP_mc = (char*)"~mnguyen/Work/bTagging/histos/ppMC_ppReco_ak3PF_BjetTrig_noIPupperCut.root";
  else fileNamePP_mc = (char*)"ppMC_ppReco_ak3PF_QCDjetTrig_noIPupperCut.root";
  if(doBjets)fileNamePbPb_mc = (char*) "/d102/mnguyen/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root";
  else fileNamePbPb_mc = (char*) "/d102/mnguyen/bTaggingOutput/ntuples/PbPbQCDMC_pt30by3_ipHICalibCentWeight_noTrig.root";

  // grab ntuples
  TFile *infPbPb_mc = new TFile(fileNamePbPb_mc);
  TFile *infPP_mc = new TFile(fileNamePP_mc);
  

  // Output file
  TFile *pbpb_Unfo;
  if (isMC) pbpb_Unfo = new TFile(Form("pbpb_Unfo_%s_MC.root",algoName[algo]),"RECREATE");
  else{
    if(doBjets)pbpb_Unfo  = new TFile(Form("pbpb_Unfo_%s_jtpt%.0f_bJets.root",algoName[algo],recoJetPtCut),"RECREATE");
    else pbpb_Unfo  = new TFile(Form("pbpb_Unfo_%s_jtpt%.0f_Inc.root",algoName[algo],recoJetPtCut),"RECREATE");
  }
  // Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];
		
  // Initialize Histograms   
	
  for (int i=0;i<=nbins_cent;i++) uhist[i] = new UnfoldingHistos(i);
	
  // Initialize reweighting functions
  
  TCut dataSelection;
  TCut dataSelectionPP;
  TCut TriggerSelectionPP;
  TCut TriggerSelectionPbPb80;

  if(doBjets)dataSelection = "weight*(abs(refparton_flavorForB)==5&&abs(jteta)<2)";
  else dataSelection = "weight*(abs(jteta)<2)";
  

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

  cout <<"MC..."<<endl;	
  
  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
 
		
  // Fill PbPb MC   
  if (!useMatrixFromFile) {
    for (Long64_t jentry2=0; jentry2<dataPbPb->tJet->GetEntries();jentry2++) {
      dataPbPb->tJet->GetEntry(jentry2);
      
      // change when we switch to centrality binning
      int cBin = 0;
            
      if ( dataPbPb->refpt  < 0. ) continue;
      if ( dataPbPb->jteta  > 2. || dataPbPb->jteta < -2. ) continue;
      if ( dataPbPb->refpt<0) dataPbPb->refpt=0;
      if (doBjets && fabs(dataPbPb->refparton_flavorForB)!=5) continue;
      if (doBjets&& dataPbPb->discr_ssvHighEff<2) continue;
      if (doBjets && dataPbPb->jtptB < recoJetPtCut) continue;
      if (!doBjets && dataPbPb->jtptA < recoJetPtCut) continue;
      //if ( dataPbPb->isTrig <1) continue;
      
      if(!doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptA>120) continue;
      if(doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptB>120) continue;

      if (!isMC||jentry2 % 2 == 1) {
	if(doBjets)uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtptB,dataPbPb->weight);
	else uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtptA,dataPbPb->weight);
      }	  
      if (jentry2 % 2 == 0) {
	uhist[cBin]-> hGen->Fill(dataPbPb->refpt,dataPbPb->weight);   
	if(doBjets)uhist[cBin]-> hMeas->Fill(dataPbPb->jtptB,dataPbPb->weight);  	 
	else uhist[cBin]-> hMeas->Fill(dataPbPb->jtptA,dataPbPb->weight);  	 
	//uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-i)),dataPbPb->weight); 
	// FIXME!!!!!!  i is supposed to be a loop over centrality !!!!
	if(doBjets)uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtptB*(1.+0.02/nbins_cent*(nbins_cent-0)),dataPbPb->weight); 
	else uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtptA*(1.+0.02/nbins_cent*(nbins_cent-0)),dataPbPb->weight); 
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
      if ( doBjets && dataPP->discr_ssvHighEff<2) continue;
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
           
    // Need to match the binning carefully, please double check whenever you change the binning
    for (int i=1;i<=hMattPbPb->GetNbinsX();i++)
      {
     	float binContent = hMattPbPb->GetBinContent(i);  
	float binError = hMattPbPb->GetBinError(i); 

	if(doBjets){
	  float tagEff =hTagEffPbPb->GetBinContent(i);
	  float tagEffErr =     hTagEffPbPb->GetBinError(i);   
	  if(tagEff>0){
	    binContent/=tagEff;
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	  }
	  else cout<<" TAGEFF = 0"<<endl;	  
	}


	uhist[0]->hMeas->SetBinContent(i+uhist[0]->hMeas->FindBin(61)-1,binContent);  
	uhist[0]->hMeas->SetBinError(i+uhist[0]->hMeas->FindBin(61)-1,binError);  
	//uhist[0]->hMeas->SetBinContent(i,hMattPbPb->GetBinContent(i));  
	//uhist[0]->hMeas->SetBinError(i,hMattPbPb->GetBinError(i));  
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
	   
    // Need to match the binning carefully, please double check whenever you change the binning
    for (int i=1;i<=hMattPP->GetNbinsX();i++)
      {
      	float binContent = hMattPP->GetBinContent(i);  
	float binError = hMattPP->GetBinError(i);  
	
	if(doBjets){
	  float tagEff =hTagEffPP->GetBinContent(i);
	  float tagEffErr =     hTagEffPP->GetBinError(i);   
	  if(tagEff>0){
	    binContent/=tagEff;
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	  }
	  else cout<<" TAGEFF = 0"<<endl;	  
	}
	
	uhist[nbins_cent]->hMeas->SetBinContent(i+uhist[nbins_cent]->hMeas->FindBin(61)-1,binContent);  
	uhist[nbins_cent]->hMeas->SetBinError(i+uhist[nbins_cent]->hMeas->FindBin(61)-1,binError);  
    
	/*
	cout<<" i "<<i<<endl;
	int newBin = i+uhist[nbins_cent]->hMeas->FindBin(61)-1;

	cout<<" newBin "<<newBin<<endl;
	cout<<"hMattPP->GetBinCenter(i) "<<hMattPP->GetBinCenter(i)<<endl;
	cout<<"uhist[nbins_cent]->hMeas->GetBinCenter(newBin) "<<uhist[nbins_cent]->hMeas->GetBinCenter(newBin)<<endl;
	*/
      }

  }


   //=========================================Response Matrix========================================================= 
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
    TH1F *hProj = (TH1F*)(TH1F*)uhist[i]->hResponse->ProjectionY(Form("hProj_cent%d",i));
    
    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {  	
	double sum=hProj->GetBinContent(y);
	cout <<y<<" "<<x<<" "<<sum<<endl;
	if (sum==0) continue;
	double ratio = uhist[i]->hMeas->GetBinContent(y)/sum;
	if (uhist[i]->hMeas->GetBinContent(y)==0) ratio = 1e-100/sum;
        uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
	uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
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
	
  TCanvas * cPbPbMeas = new TCanvas("cPbPbMeas","Measurement",1200,600);
  cPbPbMeas->Divide(2,1); 
  cPbPbMeas->cd(1);
	
  for (int i=0;i<=nbins_cent;i++) {
    cPbPb->cd(i+1)->SetLogy();   
    // Do Bin-by-bin
    TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); 
    TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",90,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
		
    // Do unfolding
    //if (isMC) uhist[i]->hMeas = (TH1F*)uhist[i]->hMatrix->ProjectionY()->Clone(Form("hMeas_cent%d",i));
    prior myPrior(uhist[i]->hMatrixFit,uhist[i]->hMeas,0);
    myPrior.unfold(uhist[i]->hMeas,1);
    TH1F *hPrior;//=(TH1F*) functionHist(fPow,uhist[i]->hMeas,Form("hPrior_cent%d",i));
//    hPrior = (TH1F*)uhist[i]->hGen->Clone("hPrior");
//    hPrior = (TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));
    hPrior=(TH1F*)hMCGen->Clone("hPrior");
    removeZero(hPrior);
    TH1F *hReweighted = (TH1F*)(TH1F*)uhist[i]->hResponse->ProjectionY(Form("hReweighted_cent%d",i));
		
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    bayesianUnfold myUnfolding(uhist[i]->hMatrixFit,myPrior.hPrior,0);
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<=40;j++)
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
    uhist[i]->hReco->SetAxisRange(0,250,"X");
    uhist[i]->hReco->Draw("");   
     
    uhist[i]->hGen->SetLineWidth(2);
    uhist[i]->hGen->SetLineColor(2);
    if(isMC)uhist[i]->hGen->Draw("hist same");
    uhist[i]->hReco->Draw("same");    
    uhist[i]->hRecoBinByBin->SetMarkerStyle(28);
    uhist[i]->hRecoBinByBin->Draw("same");    
    uhist[i]->hReco->SetAxisRange(60,240);
    TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
    hReproduced->SetMarkerColor(4);
    hReproduced->SetMarkerStyle(24);
    uhist[i]->hMeas->Draw("same");    
//    hPrior->Draw("same");
    uhist[i]->hReco->SetTitle("Baysian Unfolded");
    uhist[i]->hRecoBinByBin->SetTitle("Bin-by-bin Unfolded");
		
    TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(uhist[i]->hMeas,"Measured","pl");
    leg->AddEntry(uhist[i]->hReco,"Bayesian unfolded","pl");
    leg->AddEntry(uhist[i]->hRecoBinByBin,"Bin-by-bin unfolded","pl");
    if(isMC)leg->AddEntry(uhist[i]->hGen,"Generator level truth","l");
    leg->Draw();

    cPbPbMeas->cd(i+1)->SetLogy();   
    uhist[i]->hMeas->SetAxisRange(0,240,"X");
    uhist[i]->hMeas->Draw();
    hReproduced->Draw("same");

    TLegend *leg2 = new TLegend(0.5,0.5,0.85,0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(uhist[i]->hMeas,"Measured","pl");
    leg2->AddEntry(hReproduced,"Reproduced","pl");

    leg2->Draw();
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
         
  for (int j=2;j<20;j++) {
    hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
    hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPP[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
      hRecoIterSysPP[j]->SetTitleOffset(1.4,"Y");
      hRecoIterSysPP[j]->SetTitleOffset(1.2,"X");
      hRecoIterSysPP[j]->SetAxisRange(0,2,"Y");
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
  for (int j=2;j<20;j++) {
    hRecoIterSysPbPb[j] = rebin(uhist[0]->hRecoIterSys[j],Form("hRecoIterSysPbPb_IterSys%d",j));
    hRecoIterSysPbPb[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPbPb[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPbPb[j]->Divide(hRebinPbPb_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPbPb[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
      hRecoIterSysPbPb[j]->SetTitleOffset(1.4,"Y");
      hRecoIterSysPbPb[j]->SetTitleOffset(1.2,"X");
      hRecoIterSysPbPb[j]->SetAxisRange(0,2,"Y");
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





