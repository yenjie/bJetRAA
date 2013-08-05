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


#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"

#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldBinByBin.h"

#include "utilities.h"
#include "bayesianUnfold.h"
#include "prior.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
//==============================================================================

int Unfold2(int method =1 ,int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 1,char *spectraFileName = "pbpb_spectra_akPu3PF.root",double recoJetPtCut = 60,double trackMaxPtCut = 0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
{
#ifdef __CINT__
	gSystem->Load("libRooUnfold");
#endif
	
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

        if(method==1) {
		sprintf(chmet1,"Bayes unfo");
	} else if(method==2) {
		sprintf(chmet1,"Svd unfo ");
	} else if(method==3) {
		sprintf(chmet1,"BinByBin unfo");
	}
	
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
	
	// Define RooUnfold response
	RooUnfoldResponse *response[nbins_cent+1];
	
	// Initialize Histograms   
	
	for (int i=0;i<=nbins_cent;i++) {
		uhist[i] = new UnfoldingHistos(i);
		response[i] = new RooUnfoldResponse(uhist[i]->hResMeas,uhist[i]->hResTrue);
	}
	
	// Initialize reweighting functions
		
	TCut dataSelection;
	TCut dataSelectionPP;
        TCut TriggerSelectionPP;
        TCut TriggerSelectionPbPb80;

	if (isMC) {
	         // MC closure test, no reweighting
		dataSelection = "weight*(abs(refparton_flavorForB)==5&&abs(jteta)<2)";//Form("abs(vz)<15&&trackMax/jtpt>0.01&&abs(jteta)<2&&jtpt>%.0f&&trackMax>%f",recoJetPtCut,trackMaxPtCut);
	} else {
		dataSelection = "";
		dataSelectionPP = "";
                TriggerSelectionPP = "";//"HLT_Jet60_v1";
                TriggerSelectionPbPb80 ="HLT_HIJet80_v1";
	}
		
   	
	// Read data file
	TFile *infData;
	TFile *infPP;
	
	if (isMC) {
                infData = new TFile("/d102/mnguyen/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root");   
//   	   infData = new TFile("/d100/yjlee/hiForest/btagging/matt/unfolding/RooUnfold-1.1.1/pp/ppMC_ppReco_BjetTrig_noIPupperCut.root");
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
	
	TCanvas * cInput = new TCanvas("cInput","Input",1200,800);
	cInput->Divide(3,3);
		
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
		
	RooUnfoldResponse res(uhist[0]->hResMeas,uhist[0]->hResTrue);

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
					if (data[i]->jtpt < 60) continue;

				if (!isMC||jentry2 % 2 == 1) {
					response[cBin]->Fill(data[i]->jtpt,data[i]->refpt,data[i]->weight);
					uhist[cBin]-> hMatrix->Fill(data[i]->refpt,data[i]->jtpt,data[i]->weight);
				}	  
				if (isMC&&jentry2 % 2 == 0) {
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
	//			if(dataPP[i]->pthat<boundariesPP_pthat[i] || dataPP[i]->pthat>boundariesPP_pthat[i+1]) continue;

				int subEvt=-1;
				if ( dataPP[i]->refpt<0) continue;
				if ( dataPP[i]->jteta  > 2. || dataPP[i]->jteta < -2. ) continue;
				if ( dataPP[i]->refpt<0) dataPP[i]->refpt=0;
				if ( fabs(dataPP[i]->refparton_flavorForB)!=5) continue;
				if ( dataPP[i]->discr_ssvHighEff<2) continue;
					if (dataPP[i]->jtpt < 60) continue;
					
				if (!isMC||jentry2 % 2 == 1) {
					response[nbins_cent]->Fill(dataPP[i]->jtpt,dataPP[i]->refpt,dataPP[i]->weight);
					uhist[nbins_cent]-> hMatrix->Fill(dataPP[i]->refpt,dataPP[i]->jtpt,dataPP[i]->weight);
				}	  
				if (isMC&&jentry2 % 2 == 0) {
					uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt,dataPP[i]->weight);   
					uhist[nbins_cent]-> hMeas->Fill(dataPP[i]->jtpt,dataPP[i]->weight); 
				}           
			}
		}
	}
	

	cout <<"Response Matrix..."<<endl;
	
	TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",1200,800);
	cMatrix->Divide(3,3);

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
					//if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
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
				//uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
				//uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
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
		
//		if (!useMatrixFromFile) uhist[i]->hMatrixFit = fitMatrix(uhist[i]->hResponseNorm,recoJetPtCut);
		if (!useMatrixFromFile) uhist[i]->hMatrixFit = uhist[i]->hMatrix;
		uhist[i]->hMatrixFit->SetName(Form("hMatrixFit_cent%d",i));
	}

//	TFile *infMatt = new TFile("NewFormatV3_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
	TFile *infMatt = new TFile("/d102/mnguyen/bTaggingOutput/histos/AltBinningV5_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
	TH1F *hMatt = (TH1F*) infMatt->Get("hRawBData");
	for (int i=1;i<=hMatt->GetNbinsX();i++)
	{
	   uhist[0]->hMeas->SetBinContent(i+uhist[0]->hMeas->FindBin(61)-1,hMatt->GetBinContent(i));  
	   uhist[0]->hMeas->SetBinError(i+uhist[0]->hMeas->FindBin(61)-1,hMatt->GetBinError(i));  
	}
        divideBinWidth(hMatt);

        pbpb_Unfo->cd();
	
	cout << "==================================== TEST =====================================" << endl;
	
	cout << "==================================== UNFOLD ===================================" << endl;
	
	char chmet[100]; 
	
	// ======================= Reconstructed pp and PbPb spectra =========================================================
	TCanvas * cPbPb = new TCanvas("cPbPb","PbPb",1200,800);
	cPbPb->Divide(1,1); 
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
//		removeZero(uhist[i]->hMeas);
		
		// Do unfolding
		//if (isMC) uhist[i]->hMeas = (TH1F*)uhist[i]->hMatrix->ProjectionY()->Clone(Form("hMeas_cent%d",i));
		prior myPrior(uhist[i]->hMatrixFit,uhist[i]->hMeas,0);
		myPrior.unfold(uhist[i]->hMeas,1);
		TH1F *hPrior;//=(TH1F*) functionHist(fPow,uhist[i]->hMeas,Form("hPrior_cent%d",i));
		hPrior = (TH1F*)uhist[i]->hGen->Clone("hPrior");//(TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));
//		hPrior = (TH1F*)myPrior.hPrior->Clone("hPrior");//
//		hPrior = (TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));
//		hPrior = (TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hPrior_cent%d",i));
//		hPrior->Scale(uhist[i]->hMeas->Integral(0,1000)/uhist[i]->hRecoBinByBin->Integral(0,1000));
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
		//uhist[i]->hRecoBinByBin = (TH1F*) unfold2.Hreco();
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
				//RooUnfoldBayes unfoldToy(response[i],hToy,2);
				prior myPriorToy(hMatrixToy,hToy,0.0);
				myPriorToy.unfold(hToy,1);
				bayesianUnfold myUnfoldingToy(hMatrixToy,myPriorToy.hPrior,0.0);
				myUnfoldingToy.unfold(hToy,nBayesianIter);
				RooUnfoldBinByBin unfoldToy2(response[i],hToy);
				TH1F *hRecoTmp = (TH1F*) myUnfoldingToy.hPrior->Clone();
				TH1F *hRecoTmp2 = (TH1F*) unfoldToy2.Hreco();
				
				for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
					hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
					hTmp2[j]->Fill(hRecoTmp2->GetBinContent(j));
				}
				delete hToy;
				delete hRecoTmp;
				delete hRecoTmp2;
				delete hMatrixToy;
			}
			TF1 *f = new TF1("f","[0]*TMath::Gaus(x,[1],[2])");
			for (int j=1;j<=nbins_truth;j++)
			{
				f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
				if (hTmp[j]->GetMean()>0) {
					hTmp[j]->Fit("f","LL Q ");
					hTmp[j]->Fit("f","LL Q ");
					//	       cToy->SaveAs(Form("toy/cent-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
					//     	       cout <<j<<" "<<f->GetParameter(2)<<endl;
					uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
				}	       
				f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
				if (hTmp2[j]->GetMean()>0) {
					hTmp2[j]->Fit("f","LL Q ");
					hTmp2[j]->Fit("f","LL Q ");
					//cToy->SaveAs(Form("toy/cent2-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
					//cout <<j<<" "<<f->GetParameter(2)<<endl;
					uhist[i]->hRecoBinByBin->SetBinError(j,f->GetParameter(2));
				}	       
				delete hTmp[j];
				delete hTmp2[j];
			}
			cPbPb->cd(i+1);
		}
		//cleanup(uhist[i]->hReco);
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
		uhist[i]->hReco->SetAxisRange(90,300);
		    TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
		      hReproduced->SetMarkerColor(4);
		      hReproduced->SetMarkerStyle(24);
//		hReproduced->Draw("same");
		uhist[i]->hMeas->Draw("same");    
		    
	}	     

	pbpb_Unfo->Write();

	//***************         Calculation of RAA        ************
	
	SysData systematics;
	
	
	// ncoll uncertainty
	prepareNcollUnc(nbins_recrebin,300.);
	
	
	
	// PP histograms 

	TF1 *fPow = new TF1("fPow","[0]*x^[1]* (1-x/2760*2)^[2]");
        fPow->SetParameters(1e13,-4.5,17);
        uhist[nbins_cent]->hReco->Fit("fPow","","",50,300);
	uhist[nbins_cent]->hReco->Fit("fPow","LL","",50,300);
	uhist[nbins_cent]->hReco->Fit("fPow","LL","",50,300);
	uhist[nbins_cent]->hReco->Fit("fPow","LL","",50,300);
	uhist[nbins_cent]->hReco->Fit("fPow","LL","",50,300);
	
	
       //	uhist[nbins_cent]->hReco = functionHist(fPow,uhist[nbins_cent]->hReco,"hRecoPP");	
	
	TH1F *hRebinPP         = rebin(uhist[nbins_cent]->hReco, "hRebinPP");
	TH1F *hRebinPP_Npart   = rebin_Npart(uhist[nbins_cent]->hReco, "hRebinPP_Npart");
	TH1F *hRebinBinByBinPP = rebin(uhist[nbins_cent]->hRecoBinByBin, "hRebinBinByBinPP");
	TH1F *hRecoPP          = (TH1F*)uhist[nbins_cent]->hReco->Clone("hRecoPP");
	TH1F *hMeasPP          = (TH1F*)uhist[nbins_cent]->hMeas->Clone("hMeasPP");
	TH1F *hRebinMeasPP     = rebin(uhist[nbins_cent]->hMeas, "hRebinMeasPP");
	TH1F *hRebinGenPP      = rebin(uhist[nbins_cent]->hGen, "hRebinGenPP");
	
	
	// Scale PP histograms
	hRebinPP               ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	hRebinPP_Npart	       ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	hRebinMeasPP           ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	hRecoPP                ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	hRebinBinByBinPP       ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	hMeasPP                ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	hRebinGenPP            ->Scale(1./pplumi/CorFac[6]/0.9837/1000000);
	
	
	//***************Scale Factor for trackcut correction ************
	
	
	// correction factor for PbPb = 0.998, 0.998, 0.995, 0.991,0.987,0.977, pp = 0.966
	
	//************************     Making Canvas       ***********************
	
	
	TCanvas * cPPMCclosure = new TCanvas("cPPMCclosure","cPPMCclosure",600,450);
	
	TCanvas * cIterSys = new TCanvas("cIterSys","Iteration Systematics",1200,800);
	makeMultiPanelCanvasWithGap(cIterSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cIterSysPP = new TCanvas("cIterSysPP","Iteration Systematics for PP",800,600);
	
	TCanvas * cJECSys = new TCanvas("cJECSys","JEC Systematics",1200,800);
	makeMultiPanelCanvasWithGap(cJECSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cSmearSys = new TCanvas("cSmearSys","Smear Systematics",1200,800);
	makeMultiPanelCanvasWithGap(cSmearSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cMatrixPbPb = new TCanvas("cMatrixPbPb","Response Matrix PbPb",1400,800);
	cMatrixPbPb->Divide(3,2);
	
	TCanvas * cMatrixPP = new TCanvas("cMatrixPP","Response Matrix pp",600,450);
	
	TCanvas * cMatrixPbPbRebin = new TCanvas("cMatrixPbPbRebin"," Rebinned Response Matrix PbPb",1400,800);
	cMatrixPbPbRebin->Divide(3,2);
	
	TCanvas * cMatrixPPRebin = new TCanvas("cMatrixPPRebin","Rebinned Response Matrix pp",600,450);
	
	TCanvas * cRAANpart = new TCanvas("cRAANpart","Different Methods vs Npart",650,450);
	TCanvas * cRAANpartResult = new TCanvas("cRAANpartResult","Bayesian RAA vs Npart",650,450);
	
	TCanvas * cSpecPbPb = new TCanvas("cSpecPbPb","Jet Spectra of PbPb",1200,800);
	makeMultiPanelCanvasWithGap(cSpecPbPb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cSpecPP = new TCanvas("cSpecPP","Jet Spectra of pp",600,400);
	
	TCanvas * cRAA = new TCanvas("cRAA","RAA",1200,800);
	makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAA2 = new TCanvas("cRAA2","RAA2",1200,800);
	makeMultiPanelCanvasWithGap(cRAA2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cBayesFinal = new TCanvas("cBayesFinal"," Final Bayesian RAA akPu3PF",1200,800);
	makeMultiPanelCanvasWithGap(cBayesFinal,3,2,0.01,0.01,0.16,0.2,0.04,0.04);		
	
	TCanvas * cBayesSmrCheck = new TCanvas("cBayesSmrCheck","Bayesian and Smearing sys only RAA",1200,800);
	makeMultiPanelCanvasWithGap(cBayesSmrCheck,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cBayesSmrCheckTotsys = new TCanvas("cBayesSmrCheckTotsys","Bayesian and Smearing total sys RAA",1200,800);
	makeMultiPanelCanvasWithGap(cBayesSmrCheckTotsys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cGSVDFinal = new TCanvas("cGSVDFinal","Final GSVD RAA akPu3PF",1200,800);
	makeMultiPanelCanvasWithGap(cGSVDFinal,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAAConesBayesResult = new TCanvas("cRAAConesBayesResult"," Bayesian method in different cones in akPuPF",1200,800);
	makeMultiPanelCanvasWithGap(cRAAConesBayesResult,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAABaysAKCALO = new TCanvas("cRAABaysAKCALO","Bayesian Ak3PuPF and Ak3PuCalo",1200,800);
	makeMultiPanelCanvasWithGap(cRAABaysAKCALO,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAACones = new TCanvas("cRAACones","RAA in different cone size",1200,800);
	makeMultiPanelCanvasWithGap(cRAACones,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cSys = new TCanvas("cSys","Total Systematics",1200,800);
	makeMultiPanelCanvasWithGap(cSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAAgsvdCones = new TCanvas("cRAAgsvdCones","gsvd RAA in different cone size",1200,800);
	makeMultiPanelCanvasWithGap(cRAAgsvdCones,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAAsmearbayesCalo = new TCanvas("cRAAsmearbayesCalo","smearing and bayesian RAA in ak3PuPF and ak3PuCalo",1200,800);
	makeMultiPanelCanvasWithGap(cRAAsmearbayesCalo,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cRAABinbybin = new TCanvas("cRAABinbybin","binbybin RAA in different cone size",1200,800);
	makeMultiPanelCanvasWithGap(cRAABinbybin,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	TCanvas * cBayesSmearRatio = new TCanvas("cBayesSmearRatio"," Ratio of Bayesian and Smearing",1200,800);
	makeMultiPanelCanvasWithGap(cBayesSmearRatio,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	
	
	
	//************   MC Closure Test For PP        *************
	
	
	TLine *line = new TLine(60,1,300,1);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	
	if(isMC)
	{		
		
		cout <<"  MC Closure Test For PP   " << endl;
		
		cPPMCclosure->cd();
		
		hRebinMeasPP->Divide(hRebinGenPP);
		hRebinBinByBinPP->Divide(hRebinGenPP);
		hRebinPP->Divide(hRebinGenPP);
		
		
		hRebinPP->SetAxisRange(90,300,"X");
		hRebinPP->SetAxisRange(0,2,"Y");
		hRebinPP ->SetLineColor(kBlack);
		hRebinPP ->SetMarkerColor(kBlack);
		hRebinMeasPP->SetMarkerStyle(24);
		hRebinMeasPP->SetLineColor(4);
		hRebinMeasPP->SetMarkerColor(4);
		hRebinBinByBinPP->SetMarkerStyle(33);
		hRebinBinByBinPP->SetLineColor(kRed);
		hRebinBinByBinPP->SetMarkerColor(kRed);
		
		makeHistTitle(hRebinPP,"","Jet p_{T} (GeV/c)","Reco / Truth");
		hRebinPP->GetYaxis()->SetTitleOffset(1.5);
		hRebinPP->GetXaxis()->SetTitleOffset(1.3);
		hRebinPP->Draw("");
		hRebinBinByBinPP->Draw("same");
		hRebinMeasPP->Draw("same");
		line->Draw();
		
		TLegend *legpp = myLegend(0.52,0.65,0.85,0.9);
		legpp->AddEntry(hRebinPP,"pp Bayesian","pl");
		legpp->AddEntry(hRebinBinByBinPP,"pp Bin-by-bin","pl");
		legpp->AddEntry(hRebinMeasPP,"pp no unfolding","pl");
		legpp->Draw();
		putCMSPrel(0.2,0.83,0.06);
		drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,20);
		drawText("PYTHIA",0.6,0.4,22);
		drawText("| #eta | <2 ",0.6,0.31,22);
		
		
	}	
	
	
	
	//************   Rebinned Response Matrix pp       *************
	
	cout <<" Plotting Rebinned Response Matrix pp   " << endl;
	
	cMatrixPPRebin->cd();
	gStyle->SetPadRightMargin(0.13);	
	
	TH2F *hMatrixPPRebin = new TH2F("hMatrixPPRebin","",11,1,12,11,1,12);
	hMatrixPPRebin->GetYaxis()->SetTitle("Recojet p_{T} (GeV/c)");
	hMatrixPPRebin->GetXaxis()->SetTitle("Genjet p_{T} (GeV/c)");
	
	for(int j=1;j<=BinLabelN;j++)
	{hMatrixPPRebin->GetYaxis()->SetBinLabel(j,BinLabel[j-1]);
		hMatrixPPRebin->GetXaxis()->SetBinLabel(j,BinLabel[j-1]);}
	
	
	for (int x=1;x<=BinLabelN;x++) {
		for (int y=1;y<=BinLabelN;y++) {  	
			hMatrixPPRebin->SetBinContent(x,y,uhist[nbins_cent]->hResponseNorm->GetBinContent(x+7,y+7));
			hMatrixPPRebin->SetBinError(x,y,uhist[nbins_cent]->hResponseNorm->GetBinError(x+7,y+7));
		}
	}
	
	
	hMatrixPPRebin->GetYaxis()->SetTitleOffset(1.5);
	hMatrixPPRebin->GetXaxis()->SetTitleOffset(1.3);
	//hMatrixPPRebin->GetYaxis()->SetLabelSize(10);
	//hMatrixPPRebin->GetXaxis()->SetLabelSize(10);
	hMatrixPPRebin->Draw("textcolz");
	drawText("pp ak3PF",0.2,0.83,20);
	
	
	
	 
	 //************    Iteration Sys PP  ak3PF        *************
	 
	 cout <<" Plotting  Iteration Sys  for PP  akPu3PF   " << endl;
	 cIterSysPP->cd();
	 
	 // Iteration systematics
	 TH1F *hRecoIterSysPP[100];
	 TH1F *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, "hRebinPP_tmp");
	 TLegend *legBayesianIterPP = myLegend(0.6,0.7,0.9,0.9);
	 
	 for (int j=2;j<7;j++) {
	 hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
	 hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
	 hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
	 hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
	 if (j==2){
	 makeHistTitle(hRecoIterSysPP[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
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
	 	 
	 /*
	 
	 //************   Response Matrix pp       *************
	 
	 cout <<" Plotting Response Matrix pp   " << endl;
	 
	 cMatrixPP->cd();
	 gStyle->SetPadRightMargin(0.13);	
	 TH2F *hMatrixPP = (TH2F*)uhist[nbins_cent]->hResponseNorm->Clone("hMatrixPP");
	 makeHistTitle(hMatrixPP,"","Generator Jet p_{T} (GeV/c) ","Reconstructed Jet p_{T} (GeV/c)");
	 hMatrixPP->GetYaxis()->SetTitleOffset(1.5);
	 hMatrixPP->GetXaxis()->SetTitleOffset(1.3);
	 hMatrixPP->SetAxisRange(30,500,"X");
	 hMatrixPP->SetAxisRange(30,500,"Y");
	 hMatrixPP->Draw("colz");
	 drawText("pp ak3PF",0.2,0.83,20);
	 
	 */
	 //************     Unfolding and Measured Jet Spectra and Ratio pp      *************
	 
	 
	 cout <<" Plotting Jet Spectra of pp " << endl;
	 
	 cSpecPP->cd();
	 cSpecPP->SetLogy();
	 
	 TH1F *hRebinPP_ = rebin(uhist[nbins_cent]->hReco,"hRebinPP_" );
	 TH1F *hRebinBinByBinPP_= rebin(uhist[nbins_cent]->hRecoBinByBin,"hRebinBinByBinPP_" );
	 TH1F *hRebinMeasPP_ = rebin(uhist[nbins_cent]->hMeas,"hRebinMeasPP_" );
/*
	 divideBinWidth(hRebinPP_);
	 divideBinWidth(hRebinBinByBinPP_);
	 divideBinWidth(hRebinMeasPP_);
	*/ 
	 makeHistTitle(hRebinPP_,"","Jet p_{T} (GeV/c) ","1/N dN_{jets} / dp_{T}");
	 hRebinPP_->GetYaxis()->SetTitleOffset(1.3);
	 hRebinPP_->GetXaxis()->SetTitleOffset(1.2);
	 hRebinMeasPP_->SetMarkerStyle(24);
	 hRebinBinByBinPP_->SetMarkerStyle(34);
	 hRebinBinByBinPP_->SetMarkerColor(kRed);
	 hRebinBinByBinPP_->SetLineColor(kRed);
	 hRebinPP_->Draw();
	 hRebinMeasPP_->Draw("same");
	 hRebinBinByBinPP_->Draw("same");
	 
	 TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
	 leg->AddEntry(hRebinPP_,"Bayesian","pl");
	 leg->AddEntry(hRebinBinByBinPP_,"Bin-by-bin","pl");
	 leg->AddEntry(hRebinMeasPP_,"No unfolding","pl");
	 leg->Draw();
	 
	 
	
	
	//************   All Final Result and MC Closure Test for PbPb        *************
	
	
	TGraphErrors *smearing_vs_cent = 0;
	TGraphErrors *gsvd_vs_cent = 0;
	TGraphErrors *binbybin_vs_cent = 0;
	TGraphErrors *bayes_vs_cent = 0;
	TGraphErrors *bayes_vs_cent_allpt = 0;
	
	double xerr[6]={0};
	double yv_Bayes[6], xv_Bayes[6],valErr_Bayes[6],valErrSys_Bayes[6];
	double yv_Bayes_allpt[6], xv_Bayes_allpt[6],valErr_Bayes_allpt[6];
	double yv_BinbyBin[6], xv_BinbyBin[6],valErr_BinbyBin[6];
	double xv_smr[6], yv_smr[6],valErr_Smr[6];
	double xv_gsvd[6], yv_gsvd[6],valErr_GSVD[6];
	double SmrSys_Npart[6]={0};
	
	
//	const double npart[nbins_cent]={381.29,329.41,224.28,108.12,42.04,11.43};
	const double npart[nbins_cent]={1};
	
	// calculate the maximum from all cent bins
	for (int i=0;i<nbins_cent;i++) {
		
		TH1F *hRecoRAAIterSys[10];
		TH1F *hRebinRAA = rebin(uhist[i]->hReco, Form("hRebinRAASelected_cent%d",i));
		
		for (int j=2;j<7;j++) {
			hRecoRAAIterSys[j] = rebin(uhist[i]->hRecoIterSys[j],Form("hRecoRAA_IterSys%d_cent%d",j,i));
			hRecoRAAIterSys[j]->SetLineColor(colorCode[j-2]);
			hRecoRAAIterSys[j]->SetMarkerColor(colorCode[j-2]);
			hRecoRAAIterSys[j]->Divide(hRebinRAA);
			checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoRAAIterSys[j],0,1.1);
		}
		
		
	}
	for (int i=0;i<nbins_cent;i++) {
		cRAA->cd(nbins_cent-i);
		
		TH1F *hRebinRAA = rebin(uhist[i]->hReco, Form("hRebinRAA_cent%d",i));
		TH1F *hRebinRAA_Npart = rebin_Npart(uhist[i]->hReco, Form("hRebinRAA_Npart_cent%d",i));
		TH1F *hRebinBinByBinRAA = rebin(uhist[i]->hRecoBinByBin, Form("hRebinBinByBinRAA_cent%d",i));
		TH1F *hRebinMeasRAA = rebin(uhist[i]->hMeas, Form("hRebinMeasRAA_cent%d",i));
		TH1F *hRecoRAA         = (TH1F*)uhist[i]->hReco->Clone(Form("hRecoRAA_cent%d",i));
		TH1F *hMeasRAA         = (TH1F*)uhist[i]->hMeas->Clone(Form("hMeasRAA_cent%d",i));
		TH1F *hRecoRAAJECSys   = (TH1F*)uhist[i]->hRecoJECSys->Clone(Form("hRecoRAAJECSys_cent%d",i));
		TH1F *hRecoRAASmearSys   = (TH1F*)uhist[i]->hRecoSmearSys->Clone(Form("hRecoRAASmearSys_cent%d",i));
		TH1F *hRebinRAAJECSys = rebin(hRecoRAAJECSys, Form("hRebinRAAJECSys_cent%d",i));
		TH1F *hRebinRAASmearSys = rebin(hRecoRAASmearSys, Form("hRebinRAASmearSys_cent%d",i));
		
		TLine *l = new TLine(60,1,300,1);
		l->SetLineStyle(2);
		l->SetLineWidth(2);
		TLegend *title=0;
		if(!isMC) title = myLegend(0.18,0.7,0.48,0.8);//data
		if (isMC) title = myLegend(0.18,0.35,0.48,0.45);//MC
		title->AddEntry(hRecoRAA,Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
		title->SetTextSize(0.06);
		
		// Iteration systematics
		TH1F *hRecoRAAIterSys[10];
		TLegend *legBayesianIter = myLegend(0.6,0.6,0.9,0.9);
		legBayesianIter->SetTextSize(0.042);
		
		for (int j=2;j<7;j++) {
			
			//************    Iteration Sys   akPu3PF        *************
			
			cout <<" Plotting  Iteration Sys   akPu3PF   " << endl;
			cIterSys->cd(nbins_cent-i);
			
			hRecoRAAIterSys[j] = rebin(uhist[i]->hRecoIterSys[j],Form("hRecoRAA_IterSys%d_cent%d",j,i));
			hRecoRAAIterSys[j]->SetLineColor(colorCode[j-2]);
			hRecoRAAIterSys[j]->SetMarkerColor(colorCode[j-2]);
			hRecoRAAIterSys[j]->Divide(hRebinRAA);
			if (j==2){
				makeHistTitle(hRecoRAAIterSys[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
				hRecoRAAIterSys[j]->SetAxisRange(0,2,"Y");
				hRecoRAAIterSys[j]->Draw(); 
			} else {
				hRecoRAAIterSys[j]->Draw("same");
			}
			
			checkMaximumSys(systematics.hSysIter[i],hRecoRAAIterSys[j],0,1.1);
			legBayesianIter->AddEntry(hRecoRAAIterSys[j],Form("Iteration %d",j),"pl"); 	  
		}
		if (i==nbins_cent-1) legBayesianIter->Draw();
		title->Draw();
		l->Draw();
		
		if (i==nbins_cent-1) {
			putCMSPrel(0.2,0.83,0.06);
			drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,21);
		}	
		
		if (i==nbins_cent-2)
		{   drawText("PbPb         #sqrt{s_{NN}} = 2.76 TeV",0.2,0.83,19);
		    drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
		}
		
		if (i==nbins_cent-4)
		{	drawText("Bayesian",0.68,0.83,19);
			drawText("| #eta | < 2 ",0.28,0.83,24);
		}
		 
		 DrawPanelLabel(i);
		 
		
		if (useFixedIterativeSys) {
			// always use central systematics
			systematics.hSysIter[i] = (TH1F*) systematics.hSysIter[nbins_cent]->Clone();
		}
		drawEnvelope(systematics.hSysIter[i],"hist same");
		
		
		 //************     Unfolding and Measured Jet Spectra and Ratio PbPb      *************
		 
		 cout <<" Plotting Jet Spectra of PbPb and unfo/meas ratio " << endl;
		 
		 cSpecPbPb->cd(nbins_cent-i)->SetLogy();
		 
		 
		 TH1F *hRebinPbPb_ = rebin(uhist[i]->hReco, Form("hRebinPbPb_cent%d",i));
		 TH1F *hRebinBinByBinPbPb_= rebin(uhist[i]->hRecoBinByBin, Form("hRebinBinByBinPbPb_cent%d",i));
		 TH1F *hRebinMeasPbPb_ = rebin(uhist[i]->hMeas, Form("hRebinMeasPbPb_cent%d",i));
/*
		 divideBinWidth(hRebinPbPb_);
		 divideBinWidth(hRebinBinByBinPbPb_);
		 divideBinWidth(hRebinMeasPbPb_);
	*/	 
		 makeHistTitle(hRebinPbPb_,"","Jet p_{T} (GeV/c) ","1/N dN_{jets} / dp_{T}");
		 hRebinMeasPbPb_->SetMarkerStyle(24);
		 hRebinBinByBinPbPb_->SetMarkerStyle(34);
		 hRebinBinByBinPbPb_->SetMarkerColor(kRed);
		 hRebinBinByBinPbPb_->SetLineColor(kRed);
		 hRebinPbPb_->Draw();
		 hRebinMeasPbPb_->Draw("same");
		 hRebinBinByBinPbPb_->Draw("same");
		 title->Draw();
		 
		 if (i==nbins_cent-1) {
		 TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
		 leg->AddEntry(hRebinPbPb_,"Bayesian","pl");
		 leg->AddEntry(hRebinBinByBinPbPb_,"Bin-by-bin","pl");
		 leg->AddEntry(hRebinMeasPbPb_,"No unfolding","pl");
		 leg->Draw();
		 }
		 
		 
		
		
		cRAA->cd(nbins_cent-i); 
		makeHistTitle(hRebinRAA,"","Jet p_{T} (GeV/c)","Jet R_{AA}");
		if (!isMC) {
			
			// Scale PbPb Hisotograms  and Calculate RAA
			
			hRebinRAA            ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRebinRAA_Npart      ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRebinMeasRAA        ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRebinBinByBinRAA    ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRecoRAA             ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hMeasRAA             ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRecoRAAJECSys       ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRebinRAASmearSys     ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRebinRAAJECSys       ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRecoRAASmearSys     ->Scale(1./CorFac[i]/0.9752/1.13e9/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/TAA[i]);
			hRebinMeasRAA->Divide(hRebinMeasPP);
			hRebinRAA->Divide(hRebinPP);
			hRebinRAA_Npart->Divide(hRebinPP_Npart);
			hRebinRAAJECSys->Divide(hRebinPP);
			hRebinRAASmearSys->Divide(hRebinPP);
			hRecoRAA->Divide(hRecoPP);
			hRebinBinByBinRAA->Divide(hRebinBinByBinPP);
			hMeasRAA->Divide(hMeasPP);
			hRecoRAAJECSys->Divide(hRecoPP);
			hRecoRAASmearSys->Divide(hRecoPP);
			hRebinRAA->SetYTitle("Jet R_{AA}");      
		} else {
			
			//************   MC Closure Test For PbPb        *************
			
			
			cout <<"  MC Closure Test For PbPb   " << endl;
			
			
			TH1F *hRebinGen        = rebin(uhist[i]->hGen, Form("hRebinGen_cent%d",i));
			hRebinRAA->SetYTitle("Reco / Truth");
			hRebinMeasRAA->Divide(hRebinGen);
			hRebinRAA->Divide(hRebinGen);
			hRebinRAAJECSys->Divide(hRebinGen);
			hRebinRAASmearSys->Divide(hRebinGen);
			hRecoRAA->Divide(uhist[i]->hGen);
			hRebinBinByBinRAA->Divide(hRebinGen);
			hMeasRAA->Divide(uhist[i]->hGen);
			hRecoRAAJECSys->Divide(uhist[i]->hGen);
			
		}
		hRebinRAA->SetAxisRange(90,300,"X");
		hRebinRAA->SetAxisRange(0,2,"Y");
		hRebinMeasRAA ->SetMarkerStyle(24);
		hMeasRAA ->SetLineColor(5);
		hMeasRAA ->SetMarkerColor(5);
		hRebinMeasRAA ->SetLineColor(kBlack);
		hRebinMeasRAA ->SetMarkerColor(kBlack);
		hRecoRAA->SetMarkerStyle(24);
		hRecoRAA->SetLineColor(4);
		hRecoRAA->SetMarkerColor(4);
		hRebinBinByBinRAA->SetMarkerStyle(33);
		hRebinBinByBinRAA->SetLineColor(kRed);
		hRebinBinByBinRAA->SetMarkerColor(kRed);
		
		//************   JEC Sys  akPu3PF        *************
		
		
		hRebinRAAJECSys->Divide(hRebinRAAJECSys,hRebinRAA,1,1,"B");
		hRebinRAAJECSys->SetAxisRange(0.61,1.39,"Y");
		hRebinRAAJECSys->Draw("p");
		makeHistTitle(hRebinRAAJECSys,"","Jet p_{T} (GeV/c)","Ratio",2);
		hRebinRAAJECSys->SetAxisRange(90,300,"X");
		TF1 *fPol = new TF1("fPol","[0]");
		hRebinRAAJECSys->Fit("fPol","");
		
		
		
		
		cout <<" Plotting JEC Sys  akPu3PF  " << endl;
		
		cJECSys->cd(nbins_cent-i);
		hRebinRAAJECSys->Draw("p");
		checkMaximumSys(systematics.hSysJEC[i],functionHist(fPol,systematics.hSysJEC[i],Form("hist_sysJEC_cent%d",i)));
		l->Draw();
		
		
		//************    SmearSys   akPu3PF        *************
		
		cout <<" Plotting SmearSys akPu3PF " << endl;
		
		cSmearSys->cd(nbins_cent-i);
		
		hRebinRAASmearSys->Divide(hRebinRAASmearSys,hRebinRAA,1,1,"B");
		hRebinRAASmearSys->SetAxisRange(0.61,1.39,"Y");
		hRebinRAASmearSys->Draw("p");
		makeHistTitle(hRebinRAASmearSys,"","Jet p_{T} (GeV/c)","Ratio",2);
		hRebinRAASmearSys->SetAxisRange(90,300,"X");
		hRebinRAASmearSys->Fit("fPol","");
		hRebinRAASmearSys->Draw("p");
		checkMaximumSys(systematics.hSysSmear[i],functionHist(fPol,systematics.hSysSmear[i],Form("hist_sysSmear_cent%d",i)));
		l->Draw();
		
		
		title->Draw();
		
		/*
		 
		 //************   Response Matrix PbPb       *************
		 
		 cout <<" Plotting Response Matrix PbPb   " << endl;
		 
		 cMatrixPbPb->cd(nbins_cent-i);
		 
		 TH2F *hMatrixPbPb = (TH2F*)uhist[i]->hResponseNorm->Clone("hMatrixPbPb");
		 makeHistTitle(hMatrixPbPb,"","Generator Jet p_{T} (GeV/c) ","Reconstructed Jet p_{T} (GeV/c)");
		 hMatrixPbPb->SetAxisRange(30,500,"X");
		 hMatrixPbPb->SetAxisRange(30,500,"Y");
		 hMatrixPbPb->Draw("colz");
		 title->Draw();		
		 if (i==nbins_cent-1)
		 drawText("PbPb akPu3PF",0.2,0.83,20);
		 
		 */
		
		//************   Rebinned Response Matrix PbPb       *************
		
		cout <<" Plotting Rebinned Response Matrix PbPb   " << endl;
		
		cMatrixPbPbRebin->cd(nbins_cent-i);
		gStyle->SetPadRightMargin(0.13);	
		
		TH2F *hMatrixPbPbRebin = new TH2F("hMatrixPbPbRebin","",11,1,12,11,1,12);
		hMatrixPbPbRebin->GetYaxis()->SetTitle("Recojet p_{T} (GeV/c)");
		hMatrixPbPbRebin->GetXaxis()->SetTitle("Genjet p_{T} (GeV/c)");
		
		for(int j=1;j<=BinLabelN;j++)
		{hMatrixPbPbRebin->GetYaxis()->SetBinLabel(j,BinLabel[j-1]);
			hMatrixPbPbRebin->GetXaxis()->SetBinLabel(j,BinLabel[j-1]);
		}
		
		for (int x=1;x<=BinLabelN;x++) {
			for (int y=1;y<=BinLabelN;y++) {  	
				hMatrixPbPbRebin->SetBinContent(x,y,uhist[i]->hResponseNorm->GetBinContent(x+7,y+7));
				hMatrixPbPbRebin->SetBinError(x,y,uhist[i]->hResponseNorm->GetBinError(x+7,y+7));
			}
		}
		
		
		
		hMatrixPbPbRebin->GetYaxis()->SetTitleOffset(1.5);
		hMatrixPbPbRebin->GetXaxis()->SetTitleOffset(1.3);
		//hMatrixPbPbRebin->GetYaxis()->SetLabelSize(10);
		//hMatrixPbPbRebin->GetXaxis()->SetLabelSize(10);
		hMatrixPbPbRebin->Draw("textcolz");
		title->Draw();	
		
		if (i==nbins_cent-1)
			drawText("PbPb akPu3PF",0.2,0.83,20);
		
		
		
		
		//************    Different methods RAA in akPu3PF        *************
		
		
		cout <<" Plotting different methods RAA in akPu3PF " << endl;
		
		cRAA->cd(nbins_cent-i);
		hRebinRAA->Draw();
		if(!isMC)
		{systematics.calcTotalSys(i);
		 systematics.Draw(hRebinRAA,i);
		
		}
		
		DrawPanelLabel(i);
		
		////// Put Correlated and Uncorrelated Uncertainties in the Bayesian Unfolding histogram
		TH1F *hRebinRAA_corr = new TH1F(*hRebinRAA);
		hRebinRAA_corr->SetName(Form("hRebinRAA_corr_cent",i));
		for (Int_t j = 1; j < hRebinRAA_corr->GetNbinsX() + 1; j++) {
			Float_t x = hRebinRAA_corr->GetXaxis()->GetBinCenter(j);
			Float_t s = 0;

			switch (algo) {
			case 2:
				break;
			case 3:
				switch (i) {
				case 0:
					s = 2.044234807705282 + 0.0028322876264653307*x + 3.2827115118940787e-6*pow(x,2);
					break;
				case 1:
					s = 0.9941448362655844 + 0.011261508479420763*x - 0.000010142050572618295*pow(x,2);
					break;
				case 2:
					s = 0.6746697572381065 + 0.02023907788830483*x - 0.00004519177768958705*pow(x,2);
					break;
				case 3:
					s = 0.622822048623053 + 0.01784370943830341*x - 0.00003605323648754686*pow(x,2);
					break;
				case 4:
					s = 0.7997724959427492 + 0.013842026536881468*x - 0.000020916084265574625*pow(x,2);
					break;
				case 5:
					s = 0.5639310618959498 + 0.015922843743525164*x - 0.000030872114735102704*pow(x,2);
					break;
				}
				break;
			case 4:
				break;
			}
			hRebinRAA_corr->SetBinError(j, hRebinRAA_corr->GetBinError(j) * s);
		}
		hRebinRAA_corr->SetLineColor(kPink-4);
		hRebinRAA_corr->SetLineWidth(5);
		//hRebinRAA_corr->Draw("e0x0same");
		//hRebinRAA_corr->Draw("same");
		
		
		
		
		title->Draw();
		
		TGraphErrors *smearing, *smearing_sys;
		TGraphErrors *gsvd;
		if (!isMC) {
			smearing = new TGraphErrors(Form("data/smearing/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			smearing->SetName("smearing");
			smearing->SetMarkerStyle(21);
			smearing->SetMarkerColor(kBlue);
			smearing->SetLineColor(kBlue);
			smearing->Draw("p same");
			
			
			
			// ****************Pawan's Editting******************************
			
			smearing_sys = new TGraphErrors(Form("data/smearing/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			smearing_sys->SetName("smearing_sys");
			smearing_sys->SetMarkerStyle(21);
			smearing_sys->SetMarkerColor(kBlue);
			smearing_sys->SetLineColor(kBlue);
			

			double fx=0,fy=0;          
			for(int ibin=0;ibin<smearing->GetN();ibin++){
				smearing->GetPoint(ibin, fx, fy);
				double erSy = smearing->GetErrorX(ibin);  //! systematics 
				double erSt = smearing->GetErrorY(ibin);  //! statistics
				smearing->SetPointError(ibin,0,erSt);     //! only statistical 
				smearing_sys->SetPointError(ibin,0,erSy); //! only systematic
				
				TBox *b = new TBox(hRebinRAA_corr->GetBinLowEdge(ibin+1),fy-erSy,hRebinRAA_corr->GetBinLowEdge(ibin+2),fy+erSy);
				b->SetFillColor(kViolet+6);
				b->SetFillStyle(3006);
				b->SetLineColor(kViolet+6);
				//b->Draw();
			}
			// smearing_sys->Draw("p same");
			
			SmrSys_Npart[i] = smearing_sys->GetErrorY(1);  //! Npart systematics
			
			
			// ****************   Pawan's Editting End    ******************************
		}
		
			hRebinRAA->SetMarkerStyle(20);
			hRebinRAA->SetMarkerColor(kBlack);
			hRebinRAA->Draw("same");
			hRebinMeasRAA->Draw("same");
			
			
	if (!isMC) {
			gsvd = new TGraphErrors(Form("data/gsvd/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			gsvd->SetMarkerStyle(34);
			gsvd->SetMarkerColor(kGreen+3);
			gsvd->SetLineColor(kGreen+3);
			gsvd->Draw("p same");
		}
	
		hRebinBinByBinRAA->Draw("same");
		
		
		title->Draw();
		l->Draw();
		
		if (i==nbins_cent-1) {
			TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
			leg->SetTextSize(0.05);
			leg->AddEntry(hRebinRAA,"Bayesian","pl");
			leg->AddEntry(hRebinBinByBinRAA,"Bin-by-bin","pl");
			if (!isMC){
				leg->AddEntry(smearing,"Smearing","pl");
				leg->AddEntry(gsvd,"GSVD","pl");
			}
			leg->AddEntry(hRebinMeasRAA,"No unfolding","pl");
			leg->Draw();
			putCMSPrel(0.2,0.83,0.06);
			drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,20);
		}	
		//uhist[i]->hMeas->Write();
		
		if (i==nbins_cent-2 &&!isMC)
		{   drawText("PbPb        #sqrt{s_{NN}} = 2.76 TeV",0.2,0.83,19);
			drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
		}
		if (i==nbins_cent-2 && isMC)
			
		{   drawText("PYTHIA+HYDJET",0.5,0.83,22);
			drawText("| #eta | <2 ",0.5,0.75,22);
			
		}
		if (i==nbins_cent-4 &&!isMC)
			drawText("| #eta | < 2 ",0.28,0.83,24);
		
       
		
		
		//************     Total Sys       *************
		
		
		TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
		title_->AddEntry(hRecoRAA,Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
		title_->SetTextSize(0.06);
		
		
		cout <<" Plotting Total Systematics " << endl;
		
		cSys->cd(nbins_cent-i);
		systematics.DrawComponent(i);
		
		DrawPanelLabel(i);
		
		
		if (i==nbins_cent-1) {
			putCMSPrel(0.5,0.2,0.06);
			
		}	
		
		if (i==nbins_cent-2)
		{   drawText("PbPb  #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,19);
		    drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
		    drawText("Anti-k_{T} PF Jets   R = 0.3",0.4,0.21,19);
		}
		
		title_->Draw();
		l->Draw();
		
		
		
		
		//***********************************    If not a MC closure Test Continue on the following      ************************************************ 
		
		
		if(!isMC)
		{
			
			
			//************   Final Plot for Bayesian       *************
			
			cout <<" Plotting Final Plot RAA for Bayesian  in akPu3PF " << endl;
			
			cBayesFinal->cd(nbins_cent-i);
			
			TLine *linedum = new TLine(60,1,300,1);
			//linedum->SetLineColor(kGray);
			//For Gunther's Color Systematics Band Peference
			linedum->SetLineColor(5);
			linedum->SetLineWidth(10);
			
			TLine *linetmp = new TLine(90,1,312,1);
			linetmp->SetLineStyle(2);
			linetmp->SetLineWidth(2);
			
			TH1F *hdum0 =new TH1F("hdum0","",20,90,312);
			makeHistTitle(hdum0,"","Jet p_{T} (GeV/c)","Jet R_{AA}");
			hdum0->SetAxisRange(0,2,"Y");
			hdum0->Draw();
			
			DrawPanelLabel(i);
			
			
			hRebinRAA->Draw("same");
			systematics.calcTotalSys(i);
			systematics.Draw(hRebinRAA,i);
			hRebinRAA_corr->Draw("e0x0same");
			title->Draw();
			hRebinRAA->Draw("same");
			
			tTAAerr[i]->Draw("3");
			title->Draw();
			linetmp->Draw();
			
			if (i==nbins_cent-1) {
				drawText("Bayesian",0.68,0.83,19);
				putCMSPrel(0.2,0.83,0.06);
				drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,21);
			}	
			//uhist[i]->hMeas->Write();
			
			if (i==nbins_cent-2)
			{   drawText("PbPb         #sqrt{s_{NN}} = 2.76 TeV",0.2,0.83,19);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
			}
			
			if (i==nbins_cent-4)
			{	TLegend *leg = myLegend(0.48,0.62,0.92,0.95);
				leg->SetHeader("    Uncertainties");
				leg->SetTextSize(0.048);
				leg->AddEntry(tTAAerr[0],"TAA + Lumi","f");
				leg->AddEntry(hRebinRAA_corr,"Total statistical","l");
				leg->AddEntry(hRebinRAA,"Uncorr statistical","l");
				leg->AddEntry(linedum,"Total systematics","l");
				leg->Draw();
				drawText("| #eta | < 2 ",0.28,0.83,24);
			}
			
			
			
			
			
			
			//************    Bayesian and Smearing own systematics only AkPu3PF      *************
			
			
			cout <<" Plotting Bayesian and Smearing own systematics only AkPu3PF " << endl;
			cBayesSmrCheck->cd(nbins_cent-i);
			hRebinRAA->Draw();
			systematics.DrawUnfoErr(hRebinRAA,i);
			hRebinRAA->Draw("same");
			
			DrawPanelLabel(i);
		
			TGraphErrors *smearing_sysSmr;
			
			smearing_sysSmr = new TGraphErrors(Form("data/smearing_smrsys/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			smearing_sysSmr->SetName("smearing_sys");
			smearing_sysSmr->SetMarkerStyle(21);
			smearing_sysSmr->SetMarkerColor(kBlue);
			smearing_sysSmr->SetLineColor(kBlue);
			
			
			double fx=0,fy=0;          
			for(int i=0;i<smearing->GetN();i++){
				smearing->GetPoint(i, fx, fy);
				double erSy = smearing_sysSmr->GetErrorY(i);  //! systematics 
				
				TBox *b = new TBox(hRebinRAA_corr->GetBinLowEdge(i+1),fy-erSy,hRebinRAA_corr->GetBinLowEdge(i+2),fy+erSy);
				b->SetFillColor(kViolet+6);
				b->SetFillStyle(3001);
				b->SetLineColor(kViolet+6);
				b->Draw();
			}
			
			smearing_sysSmr->Draw("p same");
			
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-1) {
				drawText("Individual Systematic",0.56,0.83,19);
				putCMSPrel(0.2,0.83,0.05);
				drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,18);
			}	
			//uhist[i]->hMeas->Write();
			
			if (i==nbins_cent-2)
			{   drawText("PbPb  #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
			}
			
			if (i==nbins_cent-3)
			{TLegend *leg = myLegend(0.55,0.58,0.95,0.95);
				leg->SetHeader("");
				leg->SetTextSize(0.04);
				leg->AddEntry(hRebinRAA,"Bayesian ","pl");
				leg->AddEntry(smearing_sysSmr,"Smearing","pl");
				leg->Draw();
			}
			
			
			
			//************    Bayesian and Smearing total systematics AkPu3PF      *************
			
			
			cout <<" Plotting Bayesian and Smearing total systematics  AkPu3PF " << endl;
			cBayesSmrCheckTotsys->cd(nbins_cent-i);
			hRebinRAA->Draw();
			systematics.Draw(hRebinRAA,i);
			hRebinRAA->Draw("same");
			
			DrawPanelLabel(i);
		
		
			for(int i=0;i<smearing->GetN();i++){
				smearing->GetPoint(i, fx, fy);
				double erSy = smearing_sys->GetErrorY(i);  //! systematics 
				
				TBox *b = new TBox(hRebinRAA_corr->GetBinLowEdge(i+1),fy-erSy,hRebinRAA_corr->GetBinLowEdge(i+2),fy+erSy);
				b->SetFillColor(kBlue-10);
				b->SetFillStyle(3001);
				b->SetLineColor(kBlue-10);
				b->Draw();
			}
			
			smearing_sys->Draw("p same");
			
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-1) {
				drawText("Total Systematics",0.56,0.83,19);
				putCMSPrel(0.2,0.83,0.05);
				drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,18);
			}	
			//uhist[i]->hMeas->Write();
			
			if (i==nbins_cent-2)
			{   drawText("PbPb  #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt  = 150~#mub^{-1}",0.5,0.78,15);
			}
			
			if (i==nbins_cent-3)
			{TLegend *leg = myLegend(0.55,0.58,0.95,0.95);
				leg->SetHeader("");
				leg->SetTextSize(0.04);
				leg->AddEntry(hRebinRAA,"Bayesian ","pl");
				leg->AddEntry(smearing_sys,"Smearing","pl");
				leg->Draw();
			}
			
			
			//************   Final Plot for GSVD        *************
			
			
			cout <<" Plotting Final Plot RAA for GSVD in akPu3PF " << endl;
			
			cGSVDFinal->cd(nbins_cent-i);
			
			TH1F *hdum =new TH1F("hdum","",20,90,300);
			makeHistTitle(hdum,"","Jet p_{T} (GeV/c)","Jet R_{AA}");
			hdum->SetAxisRange(90,300,"X");
			hdum->SetAxisRange(0,2,"Y");
			hdum->Draw();
			
			DrawPanelLabel(i);
		
			TGraphErrors *gsvd_;
			gsvd_ = new TGraphErrors(Form("data/gsvd/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			
			gsvd_->SetMarkerStyle(20);
			gsvd_->SetMarkerColor(kBlack);
			gsvd_->SetLineColor(kBlack);
			
			systematics.calcTotalSysNoUnfolding(i);
			systematics.DrawTGraph(gsvd_,i);
			
			gsvd_->Draw("p same");
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-1) {
				drawText("GSVD Unfolding",0.57,0.83,19);
				putCMSPrel(0.2,0.83,0.05);
				drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,18);
				
			}	
			//uhist[i]->hMeas->Write();
			
			if (i==nbins_cent-2)
			{   drawText("PbPb  #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
			}
			
			
			
			
			
			//************    bayesian and no unfolding RAA in akPu3PF       *************
			
			
			
			cout <<" Plotting bayesian and no unfolding RAA in akPu3PF " << endl;
			
			cRAA2->cd(nbins_cent-i);
			hRebinRAA->Draw();
			systematics.calcTotalSys(i);
			systematics.Draw(hRebinRAA,i);
			
			//hRebinRAA_corr->Draw("e0x0same");
			hRebinRAA->Draw("same");
			hRebinMeasRAA->Draw("same");
			
			title->Draw();
			l->Draw();
			
			DrawPanelLabel(i);
		
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
				leg->AddEntry(hRebinRAA,"Bayesian","pl");
				leg->AddEntry(hRebinMeasRAA,"No unfolding","pl");
				leg->Draw();
				
				
				drawText("Anti-k_{T} Particle Flow Jets    R = 0.3",0.2,0.23,18);
				
				putCMSPrel(0.2,0.83,0.05);
			}	
			
			if (i==nbins_cent-2)
			{   drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
			}
			
			
			
			
			//************     Bayesian RAA in different cones       *************
			
			
			cout <<" Plotting Bayesian RAA in different cones " << endl;
			
			cRAACones->cd(nbins_cent-i);
			hRebinRAA->Draw();
			
			DrawPanelLabel(i);
		
			systematics.Draw(hRebinRAA,i);
			TH1F *hRebinRAAAkPu2PF;
			TH1F *hRebinRAAAkPu4PF;
			TH1F *hRebinRAAAkPu2Calo;
			TH1F *hRebinRAAAkPu3Calo;
			TH1F *hRebinRAAAkPu4Calo;
			
			if (useSpectraFromFile) {
/*
				hRebinRAAAkPu2PF = (TH1F*)akPu2PFResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu4PF = (TH1F*)akPu4PFResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu2PF->SetMarkerStyle(21);
				hRebinRAAAkPu2PF->SetLineColor(kRed);
				hRebinRAAAkPu2PF->SetMarkerColor(kRed);
				
				hRebinRAAAkPu4PF->SetMarkerStyle(22);
				hRebinRAAAkPu4PF->SetLineColor(kBlue);
				hRebinRAAAkPu4PF->SetMarkerColor(kBlue);
				
				hRebinRAAAkPu2Calo = (TH1F*)akPu2CaloResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu3Calo = (TH1F*)akPu3CaloResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu4Calo = (TH1F*)akPu4CaloResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu2Calo->SetMarkerStyle(25);
				hRebinRAAAkPu2Calo->SetLineColor(kRed);
				hRebinRAAAkPu2Calo->SetMarkerColor(kRed);
				
				hRebinRAAAkPu4Calo->SetMarkerStyle(26);
				hRebinRAAAkPu4Calo->SetLineColor(kBlue);
				hRebinRAAAkPu4Calo->SetMarkerColor(kBlue);
				
				hRebinRAAAkPu3Calo->SetMarkerStyle(24);
				hRebinRAAAkPu3Calo->SetLineColor(kBlack);
				hRebinRAAAkPu3Calo->SetMarkerColor(kBlack);
				
				hRebinRAAAkPu2PF->Draw("same");
				hRebinRAAAkPu4PF->Draw("same");
				hRebinRAAAkPu2Calo->Draw("same");
				hRebinRAAAkPu3Calo->Draw("same");
				hRebinRAAAkPu4Calo->Draw("same");
				*/
			}
			
			hRebinRAA->Draw("same");
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
				if (useSpectraFromFile) {
					leg->AddEntry(hRebinRAAAkPu2PF,"akPu2PF","pl");
					leg->AddEntry(hRebinRAA,"akPu3PF","pl");
					leg->AddEntry(hRebinRAAAkPu4PF,"akPu4PF","pl");
					leg->AddEntry(hRebinRAAAkPu2Calo,"akPu2Calo","pl");
					leg->AddEntry(hRebinRAAAkPu3Calo,"akPu3Calo","pl");
					leg->AddEntry(hRebinRAAAkPu4Calo,"akPu4Calo","pl");
				} else {
					leg->AddEntry(hRebinRAA,algoName[algo],"pl");
				}
				leg->Draw();
				
				putCMSPrel(0.2,0.83,0.05);
			}	
			
			title->Draw();
			l->Draw();
			
			
			if (i==nbins_cent-2)
			{    drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
				drawText("Bayesian",0.2,0.83,20);
			}
			
			
			//************     Bayesian RAA in different cones  Result    *************
			
			
			cout <<" Plotting Bayesian RAA in different cones Result" << endl;
			
			cRAAConesBayesResult->cd(nbins_cent-i)->RedrawAxis();
			hdum0->Draw();
			systematics.Draw(hRebinRAA,i);
			
			DrawPanelLabel(i);
			
			
			if (useSpectraFromFile) {
/*
				hRebinRAAAkPu2PF = (TH1F*)akPu2PFResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu4PF = (TH1F*)akPu4PFResult->Get(Form("hRebinRAA_cent%d",i));
				hRebinRAAAkPu2PF->SetMarkerStyle(21);
				hRebinRAAAkPu2PF->SetLineColor(kRed);
				hRebinRAAAkPu2PF->SetMarkerColor(kRed);
				
				hRebinRAAAkPu4PF->SetMarkerStyle(34);
				hRebinRAAAkPu4PF->SetLineColor(kBlue);
				hRebinRAAAkPu4PF->SetMarkerColor(kBlue);
				
				TGraphErrors *GRebinRAAAkPu4PF = HistToTgraphShift(hRebinRAAAkPu4PF,1.2);
				TGraphErrors *GRebinRAAAkPu2PF = HistToTgraphShift(hRebinRAAAkPu2PF,-1.2);
				GRebinRAAAkPu2PF->SetMarkerStyle(21);
				GRebinRAAAkPu2PF->SetLineColor(kRed);
				GRebinRAAAkPu2PF->SetMarkerColor(kRed);
				
				GRebinRAAAkPu4PF->SetMarkerStyle(34);
				GRebinRAAAkPu4PF->SetLineColor(kBlue);
				GRebinRAAAkPu4PF->SetMarkerColor(kBlue);
				
				GRebinRAAAkPu2PF->Draw("p same");
	*/
				hRebinRAA->Draw("same");
		//		GRebinRAAAkPu4PF->Draw("p same");
				
				
				
				//hRebinRAAAkPu2PF->Draw("same");
				//hRebinRAAAkPu4PF->Draw("same");
			}
			
			
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
				if (useSpectraFromFile) {
					leg->AddEntry(hRebinRAAAkPu2PF,"R = 0.2","pl");
					leg->AddEntry(hRebinRAA,"R = 0.3","pl");
					leg->AddEntry(hRebinRAAAkPu4PF,"R = 0.4","pl");
					
				}
				leg->Draw();
				
				putCMSPrel(0.2,0.83,0.06);
				drawText("Anti-k_{T} Particle Flow Jets",0.2,0.23,21);
			}	
			
			title->Draw();
			linetmp->Draw();
			
			tTAAerr[i]->Draw("3");
			
			
			if (i==nbins_cent-2)
			{   drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,19);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
				drawText("Bayesian",0.2,0.83,20);
			}
			
			if (i==nbins_cent-4)
			{TLegend *leg = myLegend(0.48,0.62,0.92,0.95);
				leg->SetHeader("    Uncertainties");
				leg->SetTextSize(0.048);
				leg->AddEntry(tTAAerr[0],"TAA + Lumi","f");
				leg->AddEntry(hRebinRAA,"Uncorr statistical","l");
				leg->AddEntry(linedum,"Total systematics","l");
				leg->Draw();
				
				drawText("| #eta | < 2 ",0.28,0.83,24);
			}
			
			
			
			
			//************     Bayesian RAA in akPu3Calo and akPu3PF       *************
			
			
			cout <<" Plotting Bayesian RAA in different cones " << endl;
			
			cRAABaysAKCALO->cd(nbins_cent-i);
			hRebinRAA->Draw();
			
			DrawPanelLabel(i);
			
			if (useSpectraFromFile) {
				
				hRebinRAAAkPu3Calo->Draw("same");
			}
			
			hRebinRAA->Draw("same");
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.62,0.75,0.95,0.9);
				if (useSpectraFromFile) {
					leg->SetTextSize(0.043);
					leg->AddEntry(hRebinRAA,"Particle Flow Jets","pl");
					leg->AddEntry(hRebinRAAAkPu3Calo,"Calorimeter Jets","pl");
				} else {
					leg->AddEntry(hRebinRAA,algoName[algo],"pl");
				}
				leg->Draw();
				
				putCMSPrel(0.2,0.83,0.06);
				drawText("Anti-k_{T} Jets   R = 0.3",0.2,0.23,21);
			}	
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-2)
			{       
			        drawText("PbPb         #sqrt{s_{NN}} = 2.76 TeV",0.2,0.83,19);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
				
			}
			
			if (i==nbins_cent-4)
			{	
				drawText("| #eta | < 2 ",0.28,0.83,24);
				drawText("Bayesian",0.55,0.83,20);
			}
			
			
			
			
			//************    Smearing and Bayesian Ratio Plot   ***    Edit By Pawan and Ying*************
			
			
			cout <<" Plotting Ratio of Bayesian and Smearing " << endl;
			
			cBayesSmearRatio->cd(nbins_cent-i);
			TH1F *hRatioSMbyBayes= new TH1F;
			
			DrawPanelLabel(i);
		
			
			hRatioSMbyBayes   = (TH1F*)hRebinRAA_corr->Clone("hRatioSMbyBayes");
			//hRatioGsvdbyBayes = (TH1F*)hRebinRAA_corr->Clone("hRatioGsvdbyBayes");
			double xv;
			double v1,e1;
			for(int i=0;i<smearing->GetN();i++){
				
				smearing->GetPoint(i,xv,v1);
				e1 = smearing->GetErrorY(i);
				double ratio = v1/hRebinRAA_corr->GetBinContent(i+1);
				double err   = ratio*sqrt(pow(e1/v1,2) + pow(hRebinRAA_corr->GetBinError(i+1)/hRebinRAA_corr->GetBinContent(i+1),2));
				hRatioSMbyBayes->SetBinContent(i+1,ratio);
				hRatioSMbyBayes->SetBinError(i+1,err);
				
				//gsvd->GetPoint(i,xv,v1);
				//e1 = gsvd->GetErrorY(i);
				//ratio = v1/hRebinRAA_corr->GetBinContent(i+1);
				//err   = ratio*sqrt(pow(e1/v1,2) + pow(hRebinRAA_corr->GetBinError(i+1)/hRebinRAA_corr->GetBinContent(i+1),2));
				//hRatioGsvdbyBayes->SetBinContent(i+1,ratio);
				//hRatioGsvdbyBayes->SetBinError(i+1,err);
				
				
			}
			
			hRatioSMbyBayes->SetYTitle("Ratio of Smearing and Bayesian");
			hRatioSMbyBayes->Draw();
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-1) {
				putCMSPrel(0.2,0.83,0.05);
				drawText("Anti-k_{T} Particle Flow Jets    R = 0.3",0.2,0.23,18);
				
			}	
			
			if (i==nbins_cent-2)
			{   drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
			}
			
			
			
			
			
			//************    gsvd in different cones ***************
			
			
			
			
			
			cout <<" Plotting gsvd in different cones" << endl;
			
			cRAAgsvdCones->cd(nbins_cent-i);
			
			hdum->Draw();
			
			DrawPanelLabel(i);
		
			TGraphErrors *gsvd_akPu2PF;
			TGraphErrors *gsvd_akPu4PF;
			
			
			if (!isMC) {
				gsvd_akPu2PF = new TGraphErrors(Form("data/gsvd_akPu2PF/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
				gsvd_akPu4PF = new TGraphErrors(Form("data/gsvd_akPu4PF/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
				
				
				gsvd->Draw("p same");
				
				gsvd_akPu2PF->SetMarkerStyle(29);
				gsvd_akPu2PF->SetLineColor(kRed);
				gsvd_akPu2PF->SetMarkerColor(kRed);
				gsvd_akPu2PF->Draw("p same");
				
				gsvd_akPu4PF->SetMarkerStyle(21);
				gsvd_akPu4PF->SetLineColor(kBlue);
				gsvd_akPu4PF->SetMarkerColor(kBlue);
				gsvd_akPu4PF->Draw("p same");
				
				//hRebinMeasRAA->Draw("same");
			}
			
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
				
				leg->AddEntry(gsvd_akPu2PF,"R = 0.2","pl");
				leg->AddEntry(gsvd,"R = 0.3","pl");
				leg->AddEntry(gsvd_akPu4PF,"R = 0.4","pl");
				//leg->AddEntry(hRebinMeasRAA,"No Unfolding","pl");
				leg->Draw();
				
				drawText("Anti-k_{T} Particle Flow Jets",0.2,0.23,18);
				
				putCMSPrel(0.2,0.83,0.05);
			}	
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-2)
			{	drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
				drawText("GSVD",0.25,0.83,20);
				
			}
			
			
			
			//************    Smearing and Bayesian for akPu3Calo ***************
			
			
			cout <<" Plotting Smearing and Bayesian for akPu3Calo" << endl;
			
			cRAAsmearbayesCalo->cd(nbins_cent-i);
			
			if (useSpectraFromFile) {
			
			makeHistTitle(hRebinRAAAkPu3Calo,"","Jet p_{T} (GeV/c)","Jet R_{AA}");
			hRebinRAAAkPu3Calo->SetAxisRange(90,300,"X");
			hRebinRAAAkPu3Calo->SetAxisRange(0,2,"Y");
			
			hRebinRAAAkPu3Calo->SetMarkerStyle(20);
			hRebinRAAAkPu3Calo->SetLineColor(kRed);
			hRebinRAAAkPu3Calo->SetMarkerColor(kRed);
			hRebinRAAAkPu3Calo->Draw("p");
			systematics.calcTotalSys(i);
			systematics.Draw(hRebinRAAAkPu3Calo,i);
			
			DrawPanelLabel(i);
		
			
			
			TGraphErrors *smearing_akPu3Calo;
			
			if (!isMC) {
				smearing_akPu3Calo = new TGraphErrors(Form("data/smearing_akPu3Calo/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
				smearing_akPu3Calo->SetMarkerStyle(21);
				smearing_akPu3Calo->SetMarkerColor(kBlue);
				smearing_akPu3Calo->SetLineColor(kBlue);
				smearing_akPu3Calo->Draw("p same");
			}
			
			hRebinMeasRAA->Draw("same");
			hRebinRAAAkPu3Calo->Draw("p same");
			
			title->Draw();
			l->Draw();
			
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
				leg->AddEntry(hRebinRAAAkPu3Calo,"Bayesian","pl");
				if (!isMC)leg->AddEntry(smearing_akPu3Calo,"Smearing","pl");
				leg->AddEntry(hRebinMeasRAA,"No unfolding","pl");
				leg->Draw();
				putCMSPrel(0.2,0.83,0.05);
				
				drawText("Anti-k_{T} Calorimeter Jets    R = 0.3",0.2,0.23,18);
			}	
			
			if (i==nbins_cent-2)
			{   drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt = 150 #mub^{-1}",0.5,0.78,15);
				
			}
			
			}
			
			
			//************    Bin By Bin in different cone sizes***************
			
			cout <<" Plotting Bin By Bin in different cone sizes" << endl;
			
			cRAABinbybin->cd(nbins_cent-i);
			makeHistTitle(hRebinBinByBinRAA,"","Jet p_{T} (GeV/c)","Jet R_{AA}");
			hRebinBinByBinRAA->SetAxisRange(90,300,"X");
			hRebinBinByBinRAA->SetAxisRange(0,2,"Y");
			hRebinBinByBinRAA->Draw();
			
			DrawPanelLabel(i);
		
			
			TH1F *hRebinBinByBinRAAAkPu2PF;
			TH1F *hRebinBinByBinRAAAkPu4PF;
			TH1F *hRebinBinByBinRAAAkPu2Calo;
			TH1F *hRebinBinByBinRAAAkPu3Calo;
			TH1F *hRebinBinByBinRAAAkPu4Calo;
			
			if (useSpectraFromFile) {
			/*
				hRebinBinByBinRAAAkPu2PF = (TH1F*)akPu2PFResult->Get(Form("hRebinBinByBinRAA_cent%d",i));
				hRebinBinByBinRAAAkPu4PF = (TH1F*)akPu4PFResult->Get(Form("hRebinBinByBinRAA_cent%d",i));
				hRebinBinByBinRAAAkPu2PF->SetMarkerStyle(21);
				hRebinBinByBinRAAAkPu2PF->SetLineColor(kRed);
				hRebinBinByBinRAAAkPu2PF->SetMarkerColor(kRed);
				
				hRebinBinByBinRAAAkPu4PF->SetMarkerStyle(22);
				hRebinBinByBinRAAAkPu4PF->SetLineColor(kBlue);
				hRebinBinByBinRAAAkPu4PF->SetMarkerColor(kBlue);
				
				hRebinBinByBinRAAAkPu2Calo = (TH1F*)akPu2CaloResult->Get(Form("hRebinBinByBinRAA_cent%d",i));
				hRebinBinByBinRAAAkPu3Calo = (TH1F*)akPu3CaloResult->Get(Form("hRebinBinByBinRAA_cent%d",i));
				hRebinBinByBinRAAAkPu4Calo = (TH1F*)akPu4CaloResult->Get(Form("hRebinBinByBinRAA_cent%d",i));
				hRebinBinByBinRAAAkPu2Calo->SetMarkerStyle(25);
				hRebinBinByBinRAAAkPu2Calo->SetLineColor(kRed);
				hRebinBinByBinRAAAkPu2Calo->SetMarkerColor(kRed);
				
				hRebinBinByBinRAAAkPu4Calo->SetMarkerStyle(26);
				hRebinBinByBinRAAAkPu4Calo->SetLineColor(kBlue);
				hRebinBinByBinRAAAkPu4Calo->SetMarkerColor(kBlue);
				
				hRebinBinByBinRAAAkPu3Calo->SetMarkerStyle(24);
				hRebinBinByBinRAAAkPu3Calo->SetLineColor(kBlack);
				hRebinBinByBinRAAAkPu3Calo->SetMarkerColor(kBlack);
				
				hRebinBinByBinRAAAkPu2PF->Draw("same");
				hRebinBinByBinRAAAkPu4PF->Draw("same");
				hRebinBinByBinRAAAkPu2Calo->Draw("same");
				hRebinBinByBinRAAAkPu3Calo->Draw("same");
				hRebinBinByBinRAAAkPu4Calo->Draw("same");
				*/
			}
			
			if (i==nbins_cent-1) {
				TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
				if (useSpectraFromFile) {
					leg->AddEntry(hRebinRAAAkPu2PF,"akPu2PF","pl");
					leg->AddEntry(hRebinRAA,"akPu3PF","pl");
					leg->AddEntry(hRebinRAAAkPu4PF,"akPu4PF","pl");
					leg->AddEntry(hRebinRAAAkPu2Calo,"akPu2Calo","pl");
					leg->AddEntry(hRebinRAAAkPu3Calo,"akPu3Calo","pl");
					leg->AddEntry(hRebinRAAAkPu4Calo,"akPu4Calo","pl");
				} else {
					leg->AddEntry(hRebinRAA,algoName[algo],"pl");
				}
				leg->Draw();
				
				putCMSPrel(0.2,0.83,0.05);
			}	
			
			title->Draw();
			l->Draw();
			
			
			if (i==nbins_cent-2)
			{    drawText("PbPb #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,15);
				drawText("#int L dt= 150 #mub^{-1}",0.5,0.78,15);
				drawText("Bin-By-Bin",0.2,0.83,20);
			}
			
			//     *************  RAA vs Npart ***********************
			
			xerr[i]=0;
			yv_Bayes[i]=0, xv_Bayes[i]=0;
			yv_Bayes_allpt[i]=0, xv_Bayes_allpt[i]=0;
			yv_BinbyBin[i]=0, xv_BinbyBin[i]=0;
			xv_smr[i]=0, yv_smr[i]=0;
			xv_gsvd[i]=0, yv_gsvd[i]=0;
			
			yv_Bayes[i] = hRebinRAA->GetBinContent(1);
			xv_Bayes[i] = hRebinRAA->GetBinCenter(1);
			yv_Bayes_allpt[i] = hRebinRAA_Npart->GetBinContent(1);
			xv_Bayes_allpt[i] = hRebinRAA_Npart->GetBinCenter(1);
			valErr_Bayes[i] = hRebinRAA_corr->GetBinError(1);
			valErrSys_Bayes[i] = systematics.hSys[i]->GetBinError(1);
			yv_BinbyBin[i] = hRebinBinByBinRAA->GetBinContent(1);
			xv_BinbyBin[i] = hRebinBinByBinRAA->GetBinCenter(1);
			valErr_BinbyBin[i] = hRebinBinByBinRAA->GetBinError(1);
			
			
			smearing->GetPoint(0, xv_smr[i], yv_smr[i]);
			valErr_Smr[i] = smearing->GetErrorY(1);
			
			gsvd->GetPoint(0, xv_gsvd[i], yv_gsvd[i]);
			valErr_GSVD[i] = gsvd->GetErrorY(1);
			
			
			//     *************  Save Output to Text ***********************
			
			//TH1F *hRebinPbPbSpec = rebin(uhist[i]->hReco, Form("hRebinPbPbSpec_cent%d",i));
			//dumpDatatoTxt(Form("%.0f-%.0f",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5), hRebinPbPbSpec, systematics.hSys[i], hRebinRAA_corr, Form("%.0f-%.0f_PbPb_JetSpec.txt",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			dumpDatatoTxt(Form("%.0f-%.0f",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5), hRebinRAA, systematics.hSys[i], hRebinRAA_corr, Form("%.0f-%.0f_Bayes_ak3.txt",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			dumpDatatoTxt(Form("%.0f-%.0f",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5), hRebinMeasRAA, systematics.hSys[i], hRebinRAA_corr, Form("%.0f-%.0f_Meas_ak3.txt",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			dumpDatatoTxt(Form("%.0f-%.0f",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5), hRebinBinByBinRAA, systematics.hSys[i], hRebinRAA_corr, Form("%.0f-%.0f_BinByBin_ak3.txt",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
			
			
			
			
			
		}//***********************************   End of  If not a MC closure Test Continue on the following      ************************************************
		
		
	}   // centrality for loop end
	
	//******  RAA vs Npart plotting ******************
	
	if (!isMC) {	
		cRAANpart->cd();
		cout <<" RAA vs Npart " << endl;
		
		
		
		TLine *lineNpart = new TLine(0.,1,400.,1);
		lineNpart->SetLineStyle(2);
		lineNpart->SetLineWidth(2);
		
		TH1F *hdum00 =new TH1F("hdum00","",20,0.,400.);
		makeHistTitle(hdum00,"","N_{part}","Jet R_{AA}");
		hdum00->SetAxisRange(0,1.8,"Y");
		hdum00->GetYaxis()->SetTitleOffset(1.3);
		hdum00->GetXaxis()->SetTitleOffset(1.3);
		hdum00->Draw();
		
		smearing_vs_cent = new TGraphErrors(nbins_cent,npart,yv_smr,xerr,valErr_Smr);
		gsvd_vs_cent= new TGraphErrors(nbins_cent,npart,yv_gsvd,xerr,valErr_GSVD);
		binbybin_vs_cent= new TGraphErrors(nbins_cent,npart,yv_BinbyBin,xerr,valErr_BinbyBin);
		bayes_vs_cent= new TGraphErrors(nbins_cent,npart,yv_Bayes,xerr,valErr_Bayes);
		
		
		for (int i=0;i<nbins_cent;i++) {
			systematics.calcTotalSys(i);
			systematics.DrawNpartSys(yv_Bayes[i],i,npart[i]);
			TBox *b_npart = new TBox(npart[i]-6.,yv_smr[i]-SmrSys_Npart[i],npart[i]+6.,yv_smr[i]+SmrSys_Npart[i]);
			b_npart->SetFillColor(kViolet+6);
			b_npart->SetFillStyle(3001);
			b_npart->SetLineColor(kViolet+6);
			// b_npart->Draw();  // Draw sys box for smearing
			
		}
		
		
		
		binbybin_vs_cent->SetMarkerStyle(33);
		binbybin_vs_cent->SetMarkerColor(kRed);
		binbybin_vs_cent->SetLineColor(kRed);
		binbybin_vs_cent->Draw("p same");
		
		gsvd_vs_cent->SetMarkerStyle(34);
		gsvd_vs_cent->SetMarkerColor(kGreen+3);
		gsvd_vs_cent->SetLineColor(kGreen+3);
		gsvd_vs_cent->Draw("p same");
		
		smearing_vs_cent->SetMarkerStyle(21);
		smearing_vs_cent->SetMarkerColor(kBlue);
		smearing_vs_cent->SetLineColor(kBlue);
		smearing_vs_cent->Draw("p same");
		
		bayes_vs_cent->SetMarkerStyle(20);
		bayes_vs_cent->SetMarkerColor(kBlack);
		bayes_vs_cent->SetLineColor(kBlack);
		bayes_vs_cent->Draw("p same");
		
		
		
		
		TLegend *leg1 = myLegend(0.6,0.65,0.95,0.9);
		leg1->SetTextSize(0.04);
		leg1->AddEntry(bayes_vs_cent,"Bayesian","pl");
		leg1->AddEntry(binbybin_vs_cent,"Bin-by-bin","pl");
		if (!isMC){
			leg1->AddEntry(smearing_vs_cent,"Smearing","pl");
			leg1->AddEntry(gsvd_vs_cent,"GSVD","pl");
		}
		
		leg1->Draw();
		drawText("Anti-k_{T} Particle Flow Jets         R = 0.3",0.2,0.23,18);	
		drawText("100 < p_{T}^{jet} (GeV/c) <110 ",0.25,0.72,16);	
		drawText("PbPb   #sqrt{s_{NN}} = 2.76 TeV",0.25,0.86,16);
		drawText("#int L dt = 150 #mub^{-1}",0.25,0.79,16);
		lineNpart->Draw();
		
		
		
		
		
		
		
		//******  RAA vs Npart Final Result plotting ******************
		
		
		cRAANpartResult->cd();
		cout <<" RAA vs Npart Final Plot " << endl;
		TLine *linedum = new TLine(60,1,300,1);
		linedum->SetLineColor(kGray);
		//For Gunther's Color Systematics Band Peference
		//linedum->SetLineColor(kSpring);
		linedum->SetLineWidth(10);
		
		TBox *b = new TBox(0,0,0,0);
		b->SetFillColor(30);
		b->SetFillStyle(3001);
		b->SetLineColor(30);
		
		hdum00->Draw();
		
		for (int i=0;i<nbins_cent;i++)
			systematics.DrawNpartSys(yv_Bayes[i],i,npart[i]);
		
		//DrawNpartTAABand();
		
		bayes_vs_cent_allpt= new TGraphErrors(nbins_cent,npart,yv_Bayes_allpt,xerr,valErr_Bayes_allpt);
		
		bayes_vs_cent->SetMarkerStyle(20);
		bayes_vs_cent->SetMarkerColor(kBlack);
		bayes_vs_cent->SetLineColor(kBlack);
		bayes_vs_cent->Draw("p same");
		
		bayes_vs_cent_allpt->SetMarkerStyle(25);
		bayes_vs_cent_allpt->SetMarkerColor(kRed);
		bayes_vs_cent_allpt->SetLineColor(kRed);
		bayes_vs_cent_allpt->Draw("p same");
		
		
		drawText("Anti-k_{T} Particle Flow Jets         R = 0.3",0.2,0.23,21);
		drawText("Bayesian",0.58,0.86,19);
		drawText("PbPb   #sqrt{s_{NN}} = 2.76 TeV",0.22,0.86,19);
		drawText("#int L dt = 150 #mub^{-1}   | #eta | < 2",0.22,0.79,19);
		putCMSPrel(0.2,0.32,0.053);
		
		TLegend *leg2= myLegend(0.57,0.63,0.89,0.83);
		leg2->SetTextSize(0.04);
		//leg2->AddEntry(b,"TAA + Lumi","f");
		leg2->AddEntry(bayes_vs_cent,"Uncorr statistical","l");
		leg2->AddEntry(linedum,"Total systematics","l");
		leg2->Draw();
		
		TLegend *leg3= myLegend(0.20,0.63,0.54,0.76);
		leg3->SetTextSize(0.04);
		leg3->AddEntry(bayes_vs_cent,"100 < p_{T}^{jet} (GeV/c) < 110","pl");
		leg3->AddEntry(bayes_vs_cent_allpt,"100 < p_{T}^{jet} (GeV/c) < 300","pl");
		leg3->Draw();
		
		lineNpart->Draw();
		
				
	}
	
	
	
	//************    Plotting ends      ***************
	
	
	
	if (isMC) {
		cRAA->Update();
		cRAA->SaveAs("MCClosureTest/MCClosureTest.gif"); 
		cRAA->SaveAs("MCClosureTest/MCClosureTest.C"); 
		cRAA->SaveAs("MCClosureTest/MCClosureTest.pdf"); 
		cIterSys->Update();
		cIterSys->SaveAs("MCClosureTest/IterSys.gif");
		cIterSys->SaveAs("MCClosureTest/IterSys.C");
		cIterSys->SaveAs("MCClosureTest/IterSys.pdf");
		cPPMCclosure->Update();
		cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.gif");
		cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.C");
		cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.pdf");
		
	} else {
	
		cRAA2->Update();
		cRAA2->SaveAs(Form("result-%s/result2.gif",algoName[algo])); 
		cRAA2->SaveAs(Form("result-%s/result2.C",algoName[algo])); 
		cRAA2->SaveAs(Form("result-%s/result2.pdf",algoName[algo]));
		
		cIterSysPP->Update();
		cIterSysPP->SaveAs(Form("result-%s/IterSysPP.gif",algoName[algo]));
		cIterSysPP->SaveAs(Form("result-%s/IterSysPP.C",algoName[algo]));
		cIterSysPP->SaveAs(Form("result-%s/IterSysPP.pdf",algoName[algo]));
		cRAA->Update();
		cRAA->SaveAs(Form("result-%s/result.gif",algoName[algo])); 
		cRAA->SaveAs(Form("result-%s/result.C",algoName[algo])); 
		cRAA->SaveAs(Form("result-%s/result.pdf",algoName[algo])); 
		cBayesFinal->Update();
		cBayesFinal->SaveAs(Form("result-%s/BayesResult.gif",algoName[algo])); 
		cBayesFinal->SaveAs(Form("result-%s/BayesResult.C",algoName[algo])); 
		cBayesFinal->SaveAs(Form("result-%s/BayesResult.pdf",algoName[algo]));
		cIterSys->Update();
		cIterSys->SaveAs(Form("result-%s/IterSys.gif",algoName[algo]));
		cIterSys->SaveAs(Form("result-%s/IterSys.C",algoName[algo]));
		cIterSys->SaveAs(Form("result-%s/IterSys.pdf",algoName[algo]));	
		cRAABaysAKCALO->Update();
		cRAABaysAKCALO->SaveAs(Form("result-%s/BayesAk3PFCaloResult.gif",algoName[algo])); 
		cRAABaysAKCALO->SaveAs(Form("result-%s/BayesAk3PFCaloResult.C",algoName[algo])); 
		cRAABaysAKCALO->SaveAs(Form("result-%s/BayesAk3PFCaloResult.pdf",algoName[algo])); 
		cRAANpartResult->Update();
		cRAANpartResult->SaveAs(Form("result-%s/result_Npart_bayes.gif",algoName[algo])); 
		cRAANpartResult->SaveAs(Form("result-%s/result_Npart_bayes.C",algoName[algo])); 
		cRAANpartResult->SaveAs(Form("result-%s/result_Npart_bayes.pdf",algoName[algo])); 
		cRAANpart->Update();
		cRAANpart->SaveAs(Form("result-%s/result_Npart.gif",algoName[algo])); 
		cRAANpart->SaveAs(Form("result-%s/result_Npart.C",algoName[algo])); 
		cRAANpart->SaveAs(Form("result-%s/result_Npart.pdf",algoName[algo])); 
		cRAAConesBayesResult->Update();
		cRAAConesBayesResult->SaveAs(Form("result-%s/BayesConesResult.gif",algoName[algo]));
		cRAAConesBayesResult->SaveAs(Form("result-%s/BayesConesResult.C",algoName[algo]));
		cRAAConesBayesResult->SaveAs(Form("result-%s/BayesConesResult.pdf",algoName[algo]));
		cSys->Update();
		cSys->SaveAs(Form("result-%s/TotalSys.gif",algoName[algo]));
		cSys->SaveAs(Form("result-%s/TotalSys.C",algoName[algo]));
		cSys->SaveAs(Form("result-%s/TotalSys.pdf",algoName[algo]));
		cMatrixPPRebin->Update();
		cMatrixPPRebin->SaveAs(Form("result-%s/BayesMatrixPPRebin.gif",algoName[algo])); 
		cMatrixPPRebin->SaveAs(Form("result-%s/BayesMatrixPPRebin.C",algoName[algo])); 
		cMatrixPPRebin->SaveAs(Form("result-%s/BayesMatrixPPRebin.pdf",algoName[algo])); 
		cMatrixPbPbRebin->Update();
		cMatrixPbPbRebin->SaveAs(Form("result-%s/BayesMatrixPbPbRebin.gif",algoName[algo])); 
		cMatrixPbPbRebin->SaveAs(Form("result-%s/BayesMatrixPbPbRebin.C",algoName[algo])); 
		cMatrixPbPbRebin->SaveAs(Form("result-%s/BayesMatrixPbPbRebin.pdf",algoName[algo])); 
		cSmearSys->Update();
		cSmearSys->SaveAs(Form("result-%s/SmearSys.gif",algoName[algo])); 
		cSmearSys->SaveAs(Form("result-%s/SmearSys.C",algoName[algo])); 
		cSmearSys->SaveAs(Form("result-%s/SmearSys.pdf",algoName[algo])); 
		cBayesSmrCheck->Update();
		cBayesSmrCheck->SaveAs(Form("result-%s/BayesSmrSysCheck.gif",algoName[algo])); 
		cBayesSmrCheck->SaveAs(Form("result-%s/BayesSmrSysCheck.C",algoName[algo])); 
		cBayesSmrCheck->SaveAs(Form("result-%s/BayesSmrSysCheck.pdf",algoName[algo])); 
		cBayesSmrCheckTotsys->Update();
		cBayesSmrCheckTotsys->SaveAs(Form("result-%s/BayesSmrTotSysCheck.gif",algoName[algo])); 
		cBayesSmrCheckTotsys->SaveAs(Form("result-%s/BayesSmrTotSysCheck.C",algoName[algo])); 
		cBayesSmrCheckTotsys->SaveAs(Form("result-%s/BayesSmrTotSysCheck.pdf",algoName[algo])); 
		cBayesSmearRatio->Update();
		cBayesSmearRatio->SaveAs(Form("result-%s/RatioBayeSmear.gif",algoName[algo])); 
		cBayesSmearRatio->SaveAs(Form("result-%s/RatioBayeSmear.C",algoName[algo])); 
		cBayesSmearRatio->SaveAs(Form("result-%s/RatioBayeSmear.pdf",algoName[algo]));
		cGSVDFinal->Update();
		cGSVDFinal->SaveAs(Form("result-%s/GSVDResult.gif",algoName[algo])); 
		cGSVDFinal->SaveAs(Form("result-%s/GSVDResult.C",algoName[algo])); 
		cGSVDFinal->SaveAs(Form("result-%s/GSVDResult.pdf",algoName[algo]));
		cJECSys->Update();
		cJECSys->SaveAs(Form("result-%s/JECSys.gif",algoName[algo]));
		cJECSys->SaveAs(Form("result-%s/JECSys.C",algoName[algo]));
		cJECSys->SaveAs(Form("result-%s/JECSys.pdf",algoName[algo]));
		
		cRAACones->Update();
		cRAACones->SaveAs(Form("result-%s/resultCones.gif",algoName[algo])); 
		cRAACones->SaveAs(Form("result-%s/resultCones.C",algoName[algo])); 
		cRAACones->SaveAs(Form("result-%s/resultCones.pdf",algoName[algo])); 
		cRAAgsvdCones->Update();
		cRAAgsvdCones->SaveAs(Form("result-%s/gsvdCones.gif",algoName[algo]));
		cRAAgsvdCones->SaveAs(Form("result-%s/gsvdCones.C",algoName[algo]));
		cRAAgsvdCones->SaveAs(Form("result-%s/gsvdCones.pdf",algoName[algo]));
		cRAAsmearbayesCalo->Update();
		cRAAsmearbayesCalo->SaveAs(Form("result-%s/smearbayesCalo.gif",algoName[algo]));
		cRAAsmearbayesCalo->SaveAs(Form("result-%s/smearbayesCalo.C",algoName[algo]));
		cRAAsmearbayesCalo->SaveAs(Form("result-%s/smearbayesCalo.pdf",algoName[algo]));
		cRAABinbybin->Update();
		cRAABinbybin->SaveAs(Form("result-%s/BinbybinCones.gif",algoName[algo]));
		cRAABinbybin->SaveAs(Form("result-%s/BinbybinCones.C",algoName[algo]));
		cRAABinbybin->SaveAs(Form("result-%s/BinbybinCones.pdf",algoName[algo]));
		cMatrixPbPb->Update();
	    cMatrixPbPb->SaveAs(Form("result-%s/BayesMatrixPbPb.gif",algoName[algo])); 
		cMatrixPbPb->SaveAs(Form("result-%s/BayesMatrixPbPb.C",algoName[algo])); 
		cMatrixPbPb->SaveAs(Form("result-%s/BayesMatrixPbPb.pdf",algoName[algo])); 
		cMatrixPP->Update();
		cMatrixPP->SaveAs(Form("result-%s/BayesMatrixPP.gif",algoName[algo])); 
		cMatrixPP->SaveAs(Form("result-%s/BayesMatrixPP.C",algoName[algo])); 
		cMatrixPP->SaveAs(Form("result-%s/BayesMatrixPP.pdf",algoName[algo]));
	
		cSpecPbPb->Update();
		cSpecPbPb->SaveAs(Form("result-%s/SpecPbPb.gif",algoName[algo])); 
		cSpecPbPb->SaveAs(Form("result-%s/SpecPbPb.C",algoName[algo])); 
		cSpecPbPb->SaveAs(Form("result-%s/SpecPbPb.pdf",algoName[algo])); 
		cSpecPP->Update();
		cSpecPP->SaveAs(Form("result-%s/SpecPP.gif",algoName[algo])); 
		cSpecPP->SaveAs(Form("result-%s/SpecPP.C",algoName[algo])); 
		cSpecPP->SaveAs(Form("result-%s/SpecPP.pdf",algoName[algo])); 
		
		
		
	}
	divideBinWidth(hCent);
	pbpb_Unfo->Write();
	//  pbpb_Unfo->Close();
}





