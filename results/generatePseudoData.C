
#include <math.h>
#include <stdio.h>

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;

using namespace RooFit;

void generatePseudoData(const char *fileName="./salida.root",
                        const char *histName="jet_pt",
                        int Nevents= 100,
                        const char *fileNameOutput="./pseudoData.root") {

    char namePdf[50]; // Array in which we will store the characters to create the pdf that will generate the pseudodata for each bin
    
    Double_t n_bini_i; // The content of each bin will be stored in this Double
    
    Double_t mean;  // The mean of the sample generated from the gaussian pdf
    
    Double_t sigma;  // The sample standard deviation of the sample generated from the gaussian pdf
    
    Int_t mean_i; // The mean rounded to the nearest integer to store as the pseudodata point
    
    Int_t first =1; //To create it for the first time
 
    // Canvas on which evertyhing will be plotted
    TCanvas* Histograma_Importado = new TCanvas("Histograma_Importado","Histograma_Importado",800,800) ;
    Histograma_Importado->Divide(1,2) ;
    Histograma_Importado->cd(1) ;
    
    
   // Import the histogram from which the pseudodata will be generated
    TFile* f_pseuDS = new TFile(fileName);  // File on which the histogram from which we will generate the pseudodata is stored
    
    TH1* h_pseuDS = (TH1*)f_pseuDS->Get(histName); // Create a local Histogram from the histogram stored on the file
    
    f_pseuDS->Draw(); // Draw the original histogram from which the pseudodata will be created
    
    // Create a copy of the original Histogram on which the pseudodata will be stored
    TH1* h_pseudoData = h_pseuDS->Clone();
    
    // Create workspace to create Gaussian functions to generate the pseudodata on each point
    RooWorkspace w("w",1);
    
    // EXAMPLE OF THE PDF GENERATED FOR BIN 5 AND THE DATA GENERATED FROM IT
    
    /*// Obtain the value of the current bin and generate
    n_bini_i = h_pseuDS->GetBinContent(20) ;
    sprintf(namePdf,"Gaussian::p(x[0,%f],mu[%f],sigma[%f])",5.0*n_bini_i,n_bini_i,sqrt(n_bini_i));
    w.factory(namePdf);
   
    // Graph the PDF generated from the bin
    Histograma_Importado->cd(2) ;
    RooPlot* xframe2 = w::x.frame(Name("xframe2"),Title("PDF to generate pseudo data for bin 3")) ;
    w::p.plotOn(xframe2);
    xframe2->Draw();
    
    // Generate the events from the PDF and graph them
    RooDataSet* obsData = w::p.generate(w::x,Nevents) ;
    Histograma_Importado->cd(3) ;
    RooPlot* xframe3 = w::x.frame(Name("xframe"),Title("Pseudo Data generated from PDF constructed from bin 3")) ;
    obsData->plotOn(xframe3);
    xframe3->Draw();
    
    // Generate all the pseudodata
   */ 
    for (int i = 1; i<h_pseuDS->GetXaxis()->GetNbins(); i++) {
        
        //Obtain the value of the current bin
        n_bini_i = h_pseuDS->GetBinContent(i) ;
    
        if(n_bini_i!=0){
 		
		if(first == 1){
			n_bini_i = h_pseuDS->GetBinContent(i) ;
	    		sprintf(namePdf,"Gaussian::p(x[0,%f],mu[%f],sigma[%f])",5.0*n_bini_i,n_bini_i,sqrt(n_bini_i));
	    		w.factory(namePdf);
			first = 0;
          	}

		else{
       		     //Define a Gaussian function with n_bini_i as its mean, and sqrt(n_bini_i) as sigma
            	   	w::x.setMax(5.0*n_bini_i) ;
            		w::mu.setVal(n_bini_i);
            		w::sigma.setVal(sqrt(n_bini_i));
            		RooDataSet* obsData = w::p.generate(w::x,Nevents) ;

         		//Calculate mean of data generated
            		mean = obsData->mean(w::x) ;
            		mean_i = TMath::Nint(mean) ;
            
            		//Calculate the standard deviation of the data created
            		sigma = obsData->sigma(w::x) ;
        
         		//Store the mean and its associated error in the bin i
          		h_pseudoData->SetBinContent(i,mean_i);
            		h_pseudoData->SetBinError(i,sigma);
        	}
        }
        
        else{ //In the event that the content of the bin is zero, it has to use a diffent method to generate the pseudodata.
            
		if(first == 1){
	    		sprintf(namePdf,"Gaussian::p(x[0,%f],mu[%f],sigma[%f])",5.0,0,sqrt(1));
	    		w.factory(namePdf);
			first = 0;
          	}

		else{
            		w::x.setMax(5.0) ;
            		w::mu.setVal(n_bini_i);
      		      	w::sigma.setVal(1.0);
     		       RooDataSet* obsData = w::p.generate(w::x,Nevents) ;
     	       
     		       //Calculate mean of data generated
     		       mean = obsData->mean(w::x) ;
     		       mean_i = TMath::Nint(mean) ;
            
     		       //Store the mean and its associated error in the bin i
     		       h_pseudoData->SetBinContent(i,mean_i) ;
     		       h_pseudoData->SetBinError(i,sqrt(mean_i)) ;
     		}       
     	   }
    }
    
    // Draw the new histogram of pseudo data
    Histograma_Importado->cd(2) ;
    h_pseudoData->SetTitle("PseudoData") ;
    h_pseudoData->SetName("PseudoData") ;
    h_pseudoData->Draw() ;
    
    // Save the new hisotgram to a .root file
    TFile* f_h_pseudoData = new TFile(fileNameOutput,"RECREATE") ;
    h_pseudoData->Write() ;
    f_h_pseudoData->Close() ;
    
}
