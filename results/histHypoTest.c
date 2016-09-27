#include <string.h>
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

/*
 The purpose of this script is to generate pseudo-data by constructing a gaussian pdf with the content of each bin. The bin content is taken as the mean of the distibution, and its square root as its standard deviation. Then a number of points are sampled from the distribution and the sample mean and the sample standard deviation are then taken as the new bin content and new bin error. This procedure is then repeated for each of the bins inside the histogram given as reference.
 
 The script takes as parameters the following (Please be aware that the values specified are dummy initializations, these parameters are redefined with whatever value you give to ithem when using the scrript):
 
    dataFileName = The name of the root file which contains the histogram of the experimental data or pseudo-data.
    dataHistName = It should be the name that the histogram has inside the .root file of the data, specified in the previous parameter
    bkgFileName = The name of the root file which contains the histogram of the total background for the channel.
    bkgHistName = It should be the name that the histogram has inside the .root file of the total background,specified in the previous parameter
    sigFileName = The name of the root file which contains the histogram of the signal tested.
    sigHistName = It should be the name that the histogram has inside the .root file of the signal, specified in the previous parameter
    calculator = type of calculator to be used (i.e: Frequentist, Hybrid, asymptotic, etc.). A detailed description of each of these can be found on the "hypothesisTest" script or in the "RooStats for searches document".
    statistic = Type of test statistic to be used in calculating the p-value.
    nToys = If the frequentist calculator is to be used, then it specifies the number of samples to be drawn to carry out the hypothesis test.
    hypoTestGraphFile = name of the file where the results of null hypothesis tests will be stored.
    newHypoTest = If it is true, it means that a new root file to store the results must be created. If it is false, it means the file specified in hypoTestGraphFile must be accesed, and the reults for the new signal point must be appended in this file. It is mainly used by the hypothesisTest shell script
 */


void histHypoTest(string dataFileName="./PsD_0.root",
                  string dataHistName="PseudoData",
                  string bkgFileName="./bkg.root",
                  string bkgHistName="DiJetMass",
                  string sigFileName="./sig1.root",
                  string sigHistName="DiJetMass",
                  int calculator = 0,
                  int statistic = 3,
                  int nToys=5000,
                  const char* hypoTestGraphFile = "hypoTestGraph.root",
                  bool newHypoTest = true) {
    
////////////////////////////////////////////////////////////////////////////////////////////////
/////              IMPORT HISTOGRAMS

    TFile* f_sig = new TFile(sigFileName.c_str());
    TFile* f_bkg = new TFile(bkgFileName.c_str());
    TFile* f_data = new TFile(dataFileName.c_str());

    TH1F* h_sig = (TH1F*)f_sig->Get(sigHistName.c_str());
    TH1F* h_bkg = (TH1F*)f_bkg->Get(bkgHistName.c_str());
    TH1F* h_data = (TH1F*)f_data->Get(dataHistName.c_str());

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
    /////////
    /////////  This step is necessary to be able to carry out the hypothesis tests using the calculators of the
    /////////  Roostats package given the histograms for data, background, and signal.
        
    HistFactory::Measurement measurement("measurement", "measurement");
    measurement.SetPOI("SigXsecOverSM"); //This Parameter is the same as /mu in profile likelihood. If set to zero it means it is the null hypothesis. If it is mu=one, it is the signal + background hypothesis. It is also used to define the exclusion limits, by setting the significance at 0.05 and solving for this parameter.

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
    
    // Add the created channel to the class measurement, from measurments we can create the workspace that will contain the ModelConfig to carry out the hypotheses tests.
    measurement.AddChannel(channel);
    
    //SIDE NOTE: Remember that modelConfig is an object in which all the important information that must be given to the Null hypothesis Calculators is stored. In other words, it is an object defined so that the Calclulators can be used more easily and uniformly.
    
    //measurement.CollectHistograms(); Por alguna razon cuando corre esta linea arroja el siguiente error: Error in <TFile::Open>: no url specified

    // The following two lines create the workspace and a ModelConfig inside it.
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
