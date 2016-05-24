// macro to fit Higgs to gg spectrum
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

void histHypoTest(const char *dataFileName="./PsD_0.root",
                  const char *dataHistName="PseudoData",
                  const char *bkgFileName="./bkg.root",
                  const char *bkgHistName="DiJetMass",
                  const char *sigFileName="./sig1.root",
                  const char *sigHistName="DiJetMass",
                  int calculator = 0,
                  int statistic = 3,
                  int nToys=5000,
                  const char* hypoTestGraphFile = "hypoTestGraph.root",
                  bool newHypoTest = true) {
    
////////////////////////////////////////////////////////////////////////////////////////////////
/////              IMPORT HISTOGRAMS

    TFile* f_sig = new TFile(sigFileName);
    TFile* f_bkg = new TFile(bkgFileName);
    TFile* f_data = new TFile(dataFileName);

    TH1F* h_sig = (TH1F*)f_sig->Get(sigHistName);
    TH1F* h_bkg = (TH1F*)f_bkg->Get(bkgHistName);
    TH1F* h_data = (TH1F*)f_data->Get(dataHistName);

    /*////////////////////////////////////////////////////////////////////////////////////////////////
    //                      DRAW HISTOGRAMS IF NEEDED
    //Declarar el canvas en el que se va a trabajar
    TCanvas* Histograma_Importado = new TCanvas("Histograma_Importado","Histograma_Importado",800,800) ;
    Histograma_Importado->Divide(2,2) ;
    Histograma_Importado->cd(1) ;

    
    //Dibujar los Histogramas
    h_sig->SetFillColor(kBlue) ;
    h_sig->SetTitle("Signal");
    h_sig->Draw();
    
    Histograma_Importado->cd(2) ;
    h_bkg->SetFillColor(kRed) ;
    h_bkg->SetTitle("Background");
    h_bkg->Draw();
    
    Histograma_Importado->cd(3) ;
    h_data->SetFillColor(kBlue) ;
    h_data->SetTitle("Data");
    h_data->Draw();

    Histograma_Importado->cd(4) ;
    THStack hs("hs","test stacked histograms");
    hs.Add(h_sig);
    hs.Add(h_bkg);
    hs.Add(h_data);
    hs.Draw("nostack");
    
    //*/
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /////////                   CREATE THE MODELCONFIG THROUGH HISTFACTORY
    /////////
    /////////  REQUIREMENTS: All the histograms must have the same number of bins
    /////////  Note that its objects have a categorization similar to their physical counterparts such as SAMPLE,
    /////////  CHANNEL, MEASURMENT
    
    cout<<"3"<<endl;
    
    HistFactory::Measurement measurement("measurement", "measurement");
    measurement.SetPOI("SigXsecOverSM"); //This Parameter is the same as /mu in profile likelihood. If set to zero it means it is the null hypothesis, and if it is one it is the signal + background. It is also used to fine the exclusion limits, by setting the significance and solving for this parameter.

    // Set the luminosity.
    // It is assumed that all histograms have been scaled by luminosity, though other forms may be used, and they are specified in the HistFactory guide.
    measurement.SetLumi(1.0);
    // Set the uncertainty around the luminosity. In this cases it is gaussian
    measurement.SetLumiRelErr(0.10);
    // For this analysis it is set constant, though it need not be set as such
    measurement.AddConstantParam("Lumi");
 
    // Declare the object of HistFactory in which the histogram of the DATA will be stored
    HistFactory::Data dataHF ;
    dataHF.SetHisto(h_data) ;

    // Declare the HistFactory object in which the SIGNAL will be stored
    HistFactory::Sample signal("signal") ;
    
    signal.AddNormFactor("SigXsecOverSM", 1, 0, 30); // This asscoiates the parameter of interest to the histogram of the signal, and sets a range from /mu = 0 (Background only hypothesis), to a signal 30 times that of the original MC signal to have a good range to find its exclusion limits
    
    signal.SetHisto(h_sig) ; // Associate the SIGNAL histogram to the object of HistFactory
    signal.ActivateStatError() ; // Include the uncertainties of the original SIGNAL histogram

    // Declare the object of HistFactory in which the histogram of the BACKGROUND will be stored
    HistFactory::Sample BKG("BKG") ;
    
    BKG.SetHisto(h_bkg) ; // Associate the BACKGROUND histogram to the object of HistFactory
    BKG.ActivateStatError() ; // Include the uncertainties of the original BACKGROUND histogram

    // Create the channel in which the objects that contain the histograms will be stored
    HistFactory::Channel channel("channel1") ;
    channel.SetData(dataHF) ; // Add the DATA to the channel
    channel.AddSample(signal) ; // Add the SIGNAL to the channel
    channel.AddSample(BKG) ; // Add the BACKGROUND to the channel
    
    // Add the created channel to the class measurement, from measurmente we can create the workspace that will contain the ModelConfig to carry out the hypotheses tests
    measurement.AddChannel(channel);
    
    //measurement.CollectHistograms(); Por alguna razon cuando corre esta linea arroja el siguiente error: Error in <TFile::Open>: no url specified

    // The following to lines create the workspace and a ModelConfig inside it.
    HistFactory::HistoToWorkspaceFactoryFast h2w(measurement) ;
    RooWorkspace* w = h2w.MakeCombinedModel(measurement) ;
    
    //Use the next line if you want to see the objects such as variables, pdfs, and additional functions contained in the created workspace
    //w->Print("t") ;

    //Exports the workspace to a new .root file to be used in the hypotheses tests
    w->writeToFile("test2.root") ;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //           TIPOS O FORMAS DE CALCULAR EL P-VALUE PARA HACER EL NULL HYPOTHESIS TEST
    //
    // type = 0 Freq calculator: Generate toys using nuisance parameter at their conditional ML estimate ( θ = θ_hat) by fitting them to the observed data. Treat constraint terms in the likelihood (e.g. systematic errors) as auxiliary measurements. introduce global observables which will be varied (tossed) for each pseudo-experiment
    //
    // type = 1 Hybrid calculator: Nuisance parameters are integrated using their pdf (the constraint term) which is interpreted as a Bayesian prior
    // type = 2 Asymptotic calculator, sólo se puede con one sided limits
    // type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
    //
    //
    //           ///////////////////////////////////////////////////////////////////////////////
    //              ESTADISTICOS QUE SE PUEDEN USAR CON LAS CALCULADORAS PARA PRUEBA DE HIPOTESIS DE DESCUBRIMIENTO
    // testStatType = 0 LEP standard simple likelihood ratio L_{s+b}(mu=1)/L_{b}(mu=0)
    //              = 1 Tevatron   ratio of profiled likelihood  L_{s+b}(mu=1,nu_hat)/L_{b}(mu=0,nu'_hat)
    //              = 2 Profile Likelihood   \lambda(mu)=L_{s+b}(mu,nu_hat)/L_{b}(mu_hat,nu'_hat)
    //              = 3 Profile Likelihood one sided (i.e. = 0 if mu_hat < 0)
   
    //USAGE:        root> StandardHypoTestDemo("fileName","workspace name","S+B modelconfig name","B model name","data set name",calculator type, test statistic type, //                             number of toys)
    //
    
    char hypoTestParam[50]; //String that will take the parameters given as imput and execute the test according to these tests
    
    //Loads the macro StandardHypoTestDemo.C so it can be used later on
    gROOT->ProcessLine(".L StandardHypoTestDemo.C");
    
    if(calculator==2 || calculator==3){
        
        //Declare the parameters pass on to StandardHypoTestDemo
        sprintf(hypoTestParam,"StandardHypoTestDemo(\"test2.root\",\"combined\",\"ModelConfig\",\"\",\"obsData\",%d,%d,%d)",calculator,statistic,newHypoTest);
        
        gROOT->ProcessLine(hypoTestParam);//Asymptotic Calculator
    }
    else if(calculator==0){
        
        //Declare the parameters pass on to StandardHypoTestDemo
        sprintf(hypoTestParam,"StandardHypoTestDemo(\"test2.root\",\"combined\",\"ModelConfig\",\"\",\"obsData\",0,%d,%d,%d,\"%s\",false,0)",statistic,newHypoTest,nToys,hypoTestGraphFile);

        gROOT->ProcessLine(hypoTestParam);//Frequentist Calculator, Despues del estadistico va el numero de Toys que usa
    }
   
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //        ESTADISTICOS QUE SE PUEDEN USAR CON LAS CALCULADORAS PARA OBTENER LOS LIMITES DE EXCLUSION
    //
    // testStatType = 0 LEP  standard simple likelihood ratio L_{s+b}(mu=1)/L_{b}(mu=0)
    //              = 1 Tevatron   ratio of profiled likelihood  L_{s+b}(mu=1,nu_hat)/L_{b}(mu=0,nu'_hat)
    //              = 2 Profile Likelihood two sided  \lambda(mu)=L_{s+b}(mu,nu_hat)/L_{b}(mu_hat,nu'_hat)
    //              = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
    //              = 4 Profile Likelihood signed ( pll = -pll if mu < mu_hat)
    //              = 5 Max Likelihood Estimate as test statistic
    //              = 6 Number of observed event as test statistic
    //
    //gROOT->ProcessLine(".L StandardHypoTestInvDemo.C");
    //gROOT->ProcessLine("StandardHypoTestInvDemo(\"test2.root\",\"combined\",\"ModelConfig\",\"\",\"asimovData\",2,2)");
    
    //*/
    
    
    
}