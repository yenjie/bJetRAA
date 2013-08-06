// Standard library
#include <math.h>
#include <iostream>
#include <fstream>

// ROOT Library
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>

// Convention 
// x: Truth
// y: Measurement

class prior {
   public:

   prior(TH2F *hResponse, TH1F *hPriorDist,double a=0.1)
   {
      hResMatrix = (TH2F*)hResponse->Clone();
      hResMatrix->SetName("hResMatrix");

      hConvMatrix = (TH2F*)hResponse->Clone();
      hConvMatrix->SetName("hConvMatrix");

      hSmearMatrix = (TH2F*)hResponse->Clone();
      hSmearMatrix->SetName("hSmearMatrix");
      
      hPrior = (TH1F*)hPriorDist->Clone();
      hPrior->SetName("hPrior");

      hUnfolded = (TH1F*)hPriorDist->Clone();
      hUnfolded->SetName("hUnfold");

      nBinsM = hResponse->GetNbinsY();
      nBinsT = hResponse->GetNbinsX();
      //cout <<"Constructor"<<endl;


      for (int m=1;m<=nBinsM;m++) {
         double sum=0;
         for (int t=1;t<=nBinsT;t++) {
	    sum+=hResMatrix->GetBinContent(t,m);
         }
         if (sum==0) continue;
	 for (int t=1;t<=nBinsT;t++) {
            hResMatrix->SetBinContent(t,m,hResMatrix->GetBinContent(t,m)/sum);
            hResMatrix->SetBinError(t,m,hResMatrix->GetBinError(t,m)/sum);
         }
      }

      for (int t=1;t<=nBinsT;t++) {
         double sum=0;
         for (int m=1;m<=nBinsM;m++) {
	    sum+=hConvMatrix->GetBinContent(t,m);
         }
         if (sum==0) continue;
         for (int m=1;m<=nBinsM;m++) {
            hConvMatrix->SetBinContent(t,m,hConvMatrix->GetBinContent(t,m)/sum);
            hConvMatrix->SetBinError(t,m,hConvMatrix->GetBinError(t,m)/sum);
         }
      }

      alpha = a;
   }
   
   ~prior(){
      delete hUnfolded;
      delete hPrior;
      delete hMeasured;
      delete hReproduced;
      delete hResMatrix;
      delete hConvMatrix;
      delete hSmearMatrix;
   };
   
   int unfold(TH1F* hM, int n); // do unfolding
   int doUnfolding();
   
    
   int nBinsM;                  // Number of bins of meaurement
   int nBinsT;                  // Number of bins of truth

   double alpha;                // Weight parameter

   TH1F* hUnfolded;             // Unfolded distribution   
   TH1F* hPrior;            	// Prior Distribution (Guess)
   TH1F* hMeasured;         	// Measured distribution
   TH1F* hReproduced;		// Reproduced

   TH2F* hResMatrix;		// Response matrix    (m->t)
   TH2F* hConvMatrix;		// Convolution matrix (t->m)
   TH2F* hSmearMatrix;          // Smearing matrix

   private:

};

// ==============================================================
//    Start Unfolding 
// ==============================================================
int prior::unfold(TH1F* hM, int n)
{
   hMeasured = (TH1F*)hM->Clone();
   hMeasured->SetName("hMeasured");
   
   hReproduced = (TH1F*)hM->Clone();
   hReproduced->SetName("hProduced");
  
   if (hMeasured->GetNbinsX() != nBinsM) {
      cout <<"The # of bins doesn't match!"<<endl;
      return 0;
   }
   
   //cout <<"Unfolding start..."<<endl;
   for (int i=0; i<n; i++) {
      int flag = doUnfolding();
      if (flag==0) {
         cout <<"Error unfolding."<<endl;
         return 0;
      }	    
     // cout <<"Unfolded iteration "<<i<<endl;
   }

   return 1;
}

// ==============================================================
//    Do Unfolding
// ==============================================================
int prior::doUnfolding()
{
   // -----------------------------------------
   // Formula:
   //    (a):   Rtm = (Rmt * Pt) / sumt' (Rmt' * Pt')  
   //    (b):   Ut = Sum_m (Rtm Mm)
   //    (c):   replace Pt by Ut
   // -----------------------------------------

   // checks
   if (hResMatrix==0) {
      cout <<"hResMatrix not found!"<<endl;
      return 0;
   }
   if (hPrior==0) {
      cout <<"hPrior not found!"<<endl;
      return 0;
   }
  

   // (a)
   //cout <<"(a)"<<endl;
   for (int m=1;m<=nBinsM;m++) {
      double sumtp = 0; 
      //cout <<m<<" "<<t<<" "<<endl;
      for (int tp=1;tp<=nBinsT;tp++) {
         sumtp += hResMatrix->GetBinContent(tp,m) * hPrior->GetBinContent(tp);
      }

      for (int t=1;t<=nBinsT;t++) {
	 double a=0;
	 if (sumtp!=0) {
	    a = hResMatrix->GetBinContent(t,m) * hPrior->GetBinContent(t) / sumtp;
	    hSmearMatrix->SetBinContent(t,m,a);
	    a = hResMatrix->GetBinContent(t,m) * hPrior->GetBinError(t) / sumtp;
	    hSmearMatrix->SetBinError(t,m,a);
	 }
      }
   }

   //cout <<"(b)"<<endl;
   // (b)
   for (int t=1;t<=nBinsT;t++) {
      double b = 0;
      double bErr = 0;
      for (int m=1;m<=nBinsM;m++) {
         b += hResMatrix->GetBinContent(t,m) * hMeasured->GetBinContent(m);
	 bErr += hResMatrix->GetBinContent(t,m) * hMeasured->GetBinError(m)*hResMatrix->GetBinContent(t,m) *	 hMeasured->GetBinError(m);
      }
      hUnfolded->SetBinContent(t,b);
      hUnfolded->SetBinError(t,sqrt(bErr));
   }

   for (int t=11;t<=nBinsT;t++) {
      double b = hUnfolded->GetBinContent(t);
      double bM = hUnfolded->GetBinContent(t-1);
      double bP = hUnfolded->GetBinContent(t+1);
      if (t==1) bM=0;
      if (t==nBinsT) bP=0;
      hUnfolded->SetBinContent(t,(1-alpha)*b+(alpha/3.)*(b+bM+bP));
      hPrior->SetBinContent(t,hUnfolded->GetBinContent(t));
      hPrior->SetBinError(t,hUnfolded->GetBinError(t));
   }

   for (int m=1;m<=nBinsM;m++) {
      double sum=0;
      for (int t=1;t<=nBinsT;t++) {
         sum += hPrior->GetBinContent(t)*hConvMatrix->GetBinContent(t,m);
	 //cout <<sum<<endl;
	 //hResMatrix->SetBinContent(t,m,hSmearMatrix->GetBinContent(t,m))
      }
      hReproduced->SetBinContent(m,sum);      
   }
   
   return 1;
}






































