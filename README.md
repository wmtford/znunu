# znunu
Text files and common code for RA2b Zinv analysers

Created from a fork of tmulholland/z-team-assembly 17 Aug. 2018.

# To create a histogram root file from Bill's front-end root macro

(in src directory) (Un)comment compiler directives in RA2bZinvAnalysis.h to select:

data/MC

skim/ntuple

TreeMaker version

The corresponding configuration of the RA2bZinvAnalysis instance is
controlled via config files RA2bZinvAnalysis.cfi (for default initialization), 
plus (filename).cfg for specialization.

Run the code with

root -l -b RA2bZinvLoadClasses.C RA2bZinvDriver.C

# To make comparisons and ratios of N-1 plots

Edit and run histoRatios.py (from python directory):

python -u histoRatios.py

# To make Z mass fits, plots, and purity calculation

Edit and run histoZmassFits.py

These load Troy's omnibus RA2b.py script to do the plotting, fitting, etc.

# Instructions for Bill's integration script

To run the integration code:

cd src/

root -l RA2bin_driver.C

# Instructions for generating new dat files (this takes time)

cd python/

./getExtrapDatFiles.py sig >&! ../datFiles/DY_signal.dat

./getExtrapDatFiles.py hdp >&! ../datFiles/DY_hdp.dat

./getExtrapDatFiles.py ldp >&! ../datFiles/DY_ldp.dat

./getDoubleRatioDatFiles.py sig >&! ../datFiles/DR_signal.dat

./getDoubleRatioDatFiles.py hdp >&! ../datFiles/DR_hdp.dat

./getDoubleRatioDatFiles.py ldp >&! ../datFiles/DR_ldp.dat

# Instructions for writing the efficiencies to the effHists.root file

cd python/

./setPhotonEffs.py

./setZllEffs.py


# code to make histograms

To go here

# syncronization stuff

To go here
