For an example of the spatial correlation calculation for Tumor CD73 expression and Monocyte CD39 expression, download CD73NatComm_Coy2021_CyCIF_TMAdata.mat from (URL pending) and run example_spatcorr.m in MATLAB. Runtime should be within a couple minutes.

corr_knn2.m is the core function for computing spatial correlations using the kNN algorithm and Pearson correlations.

subpop_spatcorr_express.m and subpop_spatcorr_express_self.m compute a given spatial correlation relation for all the TMA cores, perform exponential fitting, and evaluate statistical significance.

all_spatcorrs.m and all_spatcorrs2.m compute all the different spatial correlations (e.g. between different populations and/or expression levels) that were considered in this study. 

This code has been tested on a MATLAB v2018a on a Windows 10 OS. 

_________________________________________________________________

Similar spatial correlation analyses were repeated on CHOP TMA data, which can be found in CD73NatComm_Coy2021_CyCIF_CHOPTMAdata.mat from (URL pending).

CHOPTMA_spatcorrs.m computes all the spatial correlations on the CHOP TMA subsequently considered in this study.

CHOPTMArnaproteincycif.m plots comparisons of the average CD73 and CD39 intensities measured by CyCIF relative to protein or RNA measurements made on the same CHOP specimens.

CHOPTMA_figures.m plots the spatial correlations of CHOP TMA cores that were not from control specimens. 
