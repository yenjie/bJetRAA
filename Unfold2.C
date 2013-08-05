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

int Unfold2(int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 1,char *spectraFileName = "pbpb_spectra_akPu3PF.root",double recoJetPtCut = 60,double trackMaxPtCut = 0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
{
	int useFixedIterativeSys = 1;
	int isPyquen = 0;

        bool yinglu = 0;
	bool useSkim = 1;
	
	gStyle->SetErrorX(0.5);
	gStyle->SetPaintTextFormat("3.2f");
	gStyle->SetOptLogz(1);
	gStyle->SetPadRightMargin(0.13);	
	cout<<" ------------         Unfolding         ----------           "<<endl;
	cout<<" ==============================================================================="<<endl;
	
	int nBayesianIter = 4;
	char chmet1[100];
	
	printf("Method : %s \n",chmet1);
	
	cout << "==================================== TRAIN ====================================" << endl;
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	
	// Pthat binning
	const int nbins_pthat = 1;
	Double_t boundaries_pthat[nbins_pthat+1];
	char *fileName_pthat[nbins_pthat+1];
	Double_t xsection[nbins_pthat+1];
		
	////// MC samples
	
	boundaries_pthat[0]=30;
	fileName_pthat[0]="/d102/mnguyen/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root";   
   	xsection[0]= 1.079e-02;
	
	xsection[1] = 0;
	boundaries_pthat[1]=1000;    
	
	
	// ================ pp PtBin ======================================================================
	const int nbinsPP_pthat = 1;
	Double_t boundariesPP_pthat[nbinsPP_pthat+1];
	char *fileNamePP_pthat[nbinsPP_pthat+1];
	Double_t xsectionPP[nbinsPP_pthat+1];
	
	boundariesPP_pthat[0]=0;
	fileNamePP_pthat[0]="/d100/yjlee/hiForest/btagging/matt/unfolding/RooUnfold-1.1.1/pp/ppMC_ppReco_BjetTrig_noIPupperCut.root";   
	xsectionPP[0]= 1.079e-02;
	
	xsectionPP[1] = 0;
	boundariesPP_pthat[1]=1000;
	
	//*******************lumi number for the sample***************************
	float lumi=150.;
	float pplumi=231.;
	//*************************************************************************
	
	// Output file
	TFile *pbpb_Unfo;
	if (isMC) pbpb_Unfo = new TFile(Form("pbpb_Unfo_%s_MC.root",algoName[algo]),"RECREATE");
	else pbpb_Unfo  = new TFile(Form("pbpb_Unfo_%s_jtpt%.0f_trk%.0f.root",algoName[algo],recoJetPtCut,trackMaxPtCut),"RECREATE");
	// Histograms used by RooUnfold
	UnfoldingHistos *uhist[nbins_cent+1];
		
	// Initialize Histograms   
	
	for (int i=0;i<=nbins_cent;i++) {
		uhist[i] = new UnfoldingHistos(i);
	}
	
	// Initialize reweighting functions
		
	TCut dataSelection;
	TCut dataSelectionPP;
        TCut TriggerSelectionPP;
        TCut TriggerSelectionPbPb80;

	dataSelection = "weight*(abs(refparton_flavorForB)==5&&abs(jteta)<2)";//Form("abs(vz)<15&&trackMax/jtpt>0.01&&abs(jteta)<2&&jtpt>%.0f&&trackMax>%f",recoJetPtCut,trackMaxPtCut);
		
   	
	// Read data file
	TFile *infData;
	TFile *infPP;
	
	if (isMC) {
                infData = new TFile("/d102/mnguyen/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root");   
		cout << "This is a MC closure test"<<endl;
	} else {
                infData = new TFile("/d102/mnguyen/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root");   
		cout << "This is a data analysis"<<endl;
	}
	
	TTree *tDataJet;
        tDataJet = (TTree*)infData->Get("nt");
	
	
	if (isMC) {
   	   infPP = new TFile("/d100/yjlee/hiForest/btagging/matt/unfolding/RooUnfold-1.1.1/pp/ppMC_ppReco_BjetTrig_noIPupperCut.root");
	} else {
	   infPP = new TFile("/d100/yjlee/hiForest/btagging/matt/unfolding/RooUnfold-1.1.1/pp/ppMC_ppReco_BjetTrig_noIPupperCut.root");
	}
	
	TTree *tPPJet  = (TTree*)infPP->Get("nt");
	
	
	// Setup jet data branches, basically the jet tree branches are assigned to this object when we loop over the events
	JetData *data[nbins_pthat];   // PbPb data
	JetData *dataPP[nbins_pthat]; // pp data
	
	for (int i=0;i<nbins_pthat;i++)   data[i]   = new JetData(fileName_pthat[i],"nt","nt");	
	for (int i=0;i<nbinsPP_pthat;i++) dataPP[i] = new JetData(fileNamePP_pthat[i],"nt","nt");	
	
	TFile *fSpectra(0);	
	
	if (useSpectraFromFile||useMatrixFromFile){
		fSpectra = new TFile(spectraFileName,"read");
	}
	
	// Come back to the output file dir
	pbpb_Unfo->cd();

	// Get Jet spectra from data file
	cout <<"Reading data..."<<endl;
	
	TCanvas * cInput = new TCanvas("cInput","Input",800,400);
	cInput->Divide(2,1);
		
        cout <<"Spectra..."<<endl;	
	
	for (int i=0;i<=nbins_cent;i++){
	        cout <<nbins_cent<<endl;
		TCut centCut = Form("hiBin<%.0f&&hiBin>=%.0f",boundaries_cent[i+1],boundaries_cent[i]);
		if (useSpectraFromFile) {
			uhist[i]->hMeas = (TH1F*)fSpectra->Get(Form("hMeas_cent%d",i));
		} else {
			if (!isMC) {
				tDataJet->Project(Form("hMeas_cent%d",i),"jtpt", dataSelection&&centCut&&TriggerSelectionPbPb80);
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
	tPPJet->Project("hVzPPData","vz",dataSelectionPP&&TriggerSelectionPP);
	
        cout <<"MC..."<<endl;	
		
        TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
 
		
	// Fill PbPb MC   
	if (!useMatrixFromFile) {
		for (int i=0;i<nbins_pthat;i++) {
			if (xsection[i]==0) continue;

			cout <<"Loading pthat"<<boundaries_pthat[i]
			<<" sample, cross section = "<<xsection[i]
			<< Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[i],boundaries_pthat[i+1])
			<<"   "
			<<data[i]->tJet->GetEntries()
			<<endl;
			for (Long64_t jentry2=0; jentry2<data[i]->tJet->GetEntries();jentry2++) {
				data[i]->tJet->GetEntry(jentry2);

				int cBin = hCent->FindBin(data[i]->bin)-1;
/*
				if (cBin>=nbins_cent) continue;
				if (cBin==-1) continue;
*/
				int subEvt=-1;

				if ( data[i]->refpt  < 0. ) continue;
				if ( data[i]->jteta  > 2. || data[i]->jteta < -2. ) continue;
				if ( data[i]->refpt<0) data[i]->refpt=0;
				if ( fabs(data[i]->refparton_flavorForB)!=5) continue;
				if ( data[i]->discr_ssvHighEff<2) continue;
				if ( data[i]->jtpt < 60) continue;

				if (!isMC||jentry2 % 2 == 1) {
					uhist[cBin]-> hMatrix->Fill(data[i]->refpt,data[i]->jtpt,data[i]->weight);
				}	  
				if (jentry2 % 2 == 0) {
					uhist[cBin]-> hGen->Fill(data[i]->refpt,data[i]->weight);   
					uhist[cBin]-> hMeas->Fill(data[i]->jtpt,data[i]->weight);  	 
					uhist[cBin]-> hMeasJECSys->Fill(data[i]->jtpt*(1.+0.02/nbins_cent*(nbins_cent-i)),data[i]->weight); 
				}
			}
		}
		
		// fill pp MC
		for (int i=0;i<nbinsPP_pthat;i++) {
			if (xsectionPP[i]==0) continue;
			cout <<"Loading PP pthat"<<boundariesPP_pthat[i]
			<<" sample, cross section = "<<xsectionPP[i]
			<< Form(" pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1])<<"   "
			<<dataPP[i]->tJet->GetEntries()
			<<endl;
			for (Long64_t jentry2=0; jentry2<dataPP[i]->tJet->GetEntries();jentry2++) {
				dataPP[i]->tJet->GetEntry(jentry2);

				int subEvt=-1;
				if ( dataPP[i]->refpt<0) continue;
				if ( dataPP[i]->jteta  > 2. || dataPP[i]->jteta < -2. ) continue;
				if ( dataPP[i]->refpt<0) dataPP[i]->refpt=0;
				if ( fabs(dataPP[i]->refparton_flavorForB)!=5) continue;
				if ( dataPP[i]->discr_ssvHighEff<2) continue;
				if ( dataPP[i]->jtpt < 60) continue;
					
				if (!isMC||jentry2 % 2 == 1) {
					uhist[nbins_cent]-> hMatrix->Fill(dataPP[i]->refpt,dataPP[i]->jtpt,dataPP[i]->weight);
				}	  
				if (jentry2 % 2 == 0) {
					uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt,dataPP[i]->weight);   
					uhist[nbins_cent]-> hMeas->Fill(dataPP[i]->jtpt,dataPP[i]->weight); 
				}           
			}
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
		
		uhist[i]->hResponse->Draw("col");
		
		if (!useMatrixFromFile) uhist[i]->hMatrixFit = uhist[i]->hMatrix;
		uhist[i]->hMatrixFit->SetName(Form("hMatrixFit_cent%d",i));
	}

        if (isMC==0) {
  	   // Use measured histogram from Matt & Kurt's file
	   
	   // PbPb file:
	   TFile *infMatt = new TFile("/d102/mnguyen/bTaggingOutput/histos/AltBinningV5_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
	   TH1F *hMattPbPb = (TH1F*) infMatt->Get("hRawBData");
	   divideBinWidth(hMattPbPb);
           
	   // Need to match the binning carefully, please double check whenever you change the binning
	   for (int i=1;i<=hMattPbPb->GetNbinsX();i++)
  	   {
	      uhist[0]->hMeas->SetBinContent(i+uhist[0]->hMeas->FindBin(61)-1,hMattPbPb->GetBinContent(i));  
	      uhist[0]->hMeas->SetBinError(i+uhist[0]->hMeas->FindBin(61)-1,hMattPbPb->GetBinError(i));  
 	   }

	   // pp file:
	   // The file name needs to be updated!!!!!
	   TFile *infMattPP = new TFile("NewFormatV3_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");///d102/mnguyen/bTaggingOutput/histos/AltBinningV5_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
	   TH1F *hMattPP = (TH1F*) infMattPP->Get("hRawBData");
	   divideBinWidth(hMattPP);
	   
	   // Need to match the binning carefully, please double check whenever you change the binning
	   for (int i=1;i<=hMattPP->GetNbinsX();i++)
  	   {
	      uhist[nbins_cent]->hMeas->SetBinContent(i+uhist[nbins_cent]->hMeas->FindBin(61)-1,hMattPP->GetBinContent(i));  
	      uhist[nbins_cent]->hMeas->SetBinError(i+uhist[nbins_cent]->hMeas->FindBin(61)-1,hMattPP->GetBinError(i));  
 	   }

	}   

 
        pbpb_Unfo->cd();
	
	cout << "==================================== TEST =====================================" << endl;
	
	cout << "==================================== UNFOLD ===================================" << endl;
	
	char chmet[100]; 
	
	// ======================= Reconstructed pp and PbPb spectra =========================================================
	TCanvas * cPbPb = new TCanvas("cPbPb","PbPb",1200,600);
	cPbPb->Divide(2,1); 
	cPbPb->cd(1);
	
	
	for (int i=0;i<=nbins_cent;i++) {
		cPbPb->cd(i+1)->SetLogy();   
		// Do Bin-by-bin
		TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); 
		TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
		hBinByBinCorRaw->Divide(hMCGen);
		TF1 *f = new TF1("f","[0]+[1]*x");
		hBinByBinCorRaw->Fit("f","LL ","",90,300);
		TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
		delete hBinByBinCorRaw,hMCGen;
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
			TF1 *f = new TF1("f","[0]*TMath::Gaus(x,[1],[2])");
			for (int j=1;j<=nbins_truth;j++)
			{
				f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
				if (hTmp[j]->GetMean()>0) {
					hTmp[j]->Fit("f","LL Q ");
					hTmp[j]->Fit("f","LL Q ");
					uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
				}	       
				f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
				if (hTmp2[j]->GetMean()>0) {
					hTmp2[j]->Fit("f","LL Q ");
					hTmp2[j]->Fit("f","LL Q ");
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
		uhist[i]->hReco->Draw("");    
		uhist[i]->hGen->SetLineWidth(2);
		uhist[i]->hGen->SetLineColor(2);
		uhist[i]->hGen->Draw("hist same");
		uhist[i]->hReco->Draw("same");    
		uhist[i]->hRecoBinByBin->SetMarkerStyle(28);
		uhist[i]->hRecoBinByBin->Draw("same");    
		uhist[i]->hReco->SetAxisRange(60,300);
                TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
		hReproduced->SetMarkerColor(4);
		hReproduced->SetMarkerStyle(24);
		uhist[i]->hMeas->Draw("same");    
		
		TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(uhist[i]->hMeas,"Measured","pl");
		leg->AddEntry(uhist[i]->hReco,"Bayesian unfolded","pl");
		leg->AddEntry(uhist[i]->hRecoBinByBin,"Bin-by-bin unfolded","pl");
		leg->AddEntry(uhist[i]->hGen,"Generator level truth","l");
		leg->Draw();
	}	     

	pbpb_Unfo->Write();

        SysData systematics;
         TLine *line = new TLine(60,1,250,1);

         // Iteration systematics
	 TCanvas *cIterSys = new TCanvas("cIterSys","",1200,600);
	 cIterSys->Divide(2,1);
	 cIterSys->cd(2);
         TH1F *hRecoIterSysPP[100];
         TH1F *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, "hRebinPP_tmp");
         TLegend *legBayesianIterPP = myLegend(0.4,0.7,0.9,0.9);
         legBayesianIterPP->AddEntry("","PP","");
         
         for (int j=2;j<7;j++) {
            hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
            hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
            hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
            hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
            if (j==2){
               makeHistTitle(hRecoIterSysPP[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
 	       hRecoIterSysPP[j]->SetTitleOffset(1.5,"Y");
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
         drawEnvelope(systematics.hSysIter[nbins_cent],"hist same");


	 cIterSys->cd(1);
         TH1F *hRecoIterSysPbPb[100];
         TH1F *hRebinPbPb_tmp         = rebin(uhist[0]->hReco, "hRebinPbPb_tmp");
         TLegend *legBayesianIterPbPb = myLegend(0.4,0.7,0.9,0.9);
         legBayesianIterPbPb->AddEntry("","PbPb","");
         for (int j=2;j<7;j++) {
            hRecoIterSysPbPb[j] = rebin(uhist[0]->hRecoIterSys[j],Form("hRecoIterSysPbPb_IterSys%d",j));
            hRecoIterSysPbPb[j]->SetLineColor(colorCode[j-2]);
            hRecoIterSysPbPb[j]->SetMarkerColor(colorCode[j-2]);
            hRecoIterSysPbPb[j]->Divide(hRebinPbPb_tmp);
            if (j==2){
               makeHistTitle(hRecoIterSysPbPb[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
               hRecoIterSysPbPb[j]->SetTitleOffset(1.5,"Y");
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
         drawEnvelope(systematics.hSysIter[0],"hist same");
}





