
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

void graphFun() {
    
 
    // Canvas on which evertyhing will be plotted
    TCanvas* Histograma_Importado = new TCanvas("Histograma_Importado","Histograma_Importado",800,800) ;
    Histograma_Importado->Divide(1) ;
    Histograma_Importado->cd(1) ;
    
    // Create workspace to create Gaussian functions to generate the pseudodata on each point
    RooWorkspace w("w",1);
    
    w.factory("Gaussian::p(x[-2,6],mu[2],sigma[1])");
   
    //Graph the PDF
    RooPlot* xframe2 = w::x.frame(Name("xframe2"),Title("PDF to generate pseudo data for bin 3")) ;
   // w::p.plotOn(xframe2);
    
    
    w.factory("EXPR::kpos('0.5*(1/CosH((Pi()/2)*x))',x)");
    w::kpos.plotOn(xframe2);
    
    xframe2->Draw();
}
