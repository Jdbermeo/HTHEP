# HTHEP
Scripts for Hypothesis Tests using RooStats.

These are the steps to execute the current version of script written to execute a hypothesis test wit a number n histograms for the simulated signal and m different histograms for the simulated background, and/or an experimental data histogram as input. The output of the script is a .root document containing a tree with four results from the null hypothesis tests (p_0,CL_b=p_0,CL_{s+b}=p_1, CL_s=P_1/(1-p_0). The name of the output file is hypoTestDisc.root. 

# The steps are the following

• Clone from this repository.

• Copy the .root files containing the histograms of the background to the directory named background. Make sure all histograms have the same binning and range.

• Copy the .root files containing the histograms for each simulated signal to the directory named signal. Make sure all histograms have the same binning and range as the background ones and between them.

• If an experimental data histogram is available, copy it to the data directory, if not leave it empty.

• Open the terminal and go to the directory HTHEP. Run the following command "vi hypothesisTest".

• Modify the parameters used in the script. The parameters are:

  – n events pd: If using pseudodata, declare the number of points generated from a gaussian pdf constructed around the content of the bin. If it is too high, then there will be no difference with the histogram of signal plus background.

  – If using data, declare the name of the .root file containing the data.
  
  – Specify the name that the histograms have in the .root file. If you do not know it, run the command ls() in root to find out the name of the histograms. The shell script assumes all background histograms have the same name as objects. The same applies to the entire set of signal files. Note that the name of the file is not the same as the histograms name within the .root file. They are often different.

  – calculator: specify a number according to the type of calculator you wish to use. The Hybrid calculator is not yet available.

  – ntoys: Declare the number of toys you wish to use if using the frequentist
calculator

  – Exit the shell script
  
• Run the command ”chmod 755 hypothesisTest”. Now run ”./hypothesisTest

• First of all the shell script will use the command hadd to sum all the background histograms into a histogram.

• If there is no data file, then it will create a new set of histograms by summing each of the signals with the total background.

• It will run the script generatePseudoData.C for each of the newly created histograms and generate a set of pseudodata, one file for each signal point.

• The script now will run the script histHypoTest.C. This script does the following:

  – The script takes as input the names of the files that contain the total background, the signal point being evaluated, and the pseudodata. It takes as parameter the name that each histogram takes within its respective file.

  – It creates the ModelConfig using HistFactory and taking into account the errors of each histogram.

  – It runs a modified version of Lorenzo Moretta’s script standardHypoTestDemo.C. The shell script has also passed as parameters to histHypoTest.C the required parameters to run standardHypoTestDemo.C.

  – The modified version of standardHypoTestDemo.C. saves the results of the hypothesis test to the file ”hypoTestDisc.root”. It stores the calculated values for the the p-values p0, p1, and CLS.
  
• The script deletes all unnecessary files that were created

• Access the hypoTestDisc.root file and graph the data for all the signal points from the file. One graph for the p-values of discovery, and one graph for CLS.
