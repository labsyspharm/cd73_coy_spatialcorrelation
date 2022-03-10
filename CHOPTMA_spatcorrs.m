load CD73NatComm_Coy2021_CHOPTMAdata.mat
allcores = unique(unique_core_num);%
%% tumor cd73-express vs myeloid cd39-express
popAid = tumor_id;
popBid = myeloid_id;
magA = cleandata.CD73;
magB = cleandata.CD39;
tic
[spatcorr_cd73_cd39,radfunc_cd73_cd39,parms_cd73_cd39,h_cd73_cd39] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc
%% tumor cd73-express vs all HIF1-express
popAid = tumor_id; %tumor + ambig
popBid = tumor_id==tumor_id; % all cells
magA = cleandata.CD73;
magB = cleandata.HIF1a;
tic
[spatcorr_cd73_hif1,radfunc_cd73_hif1,parms_cd73_hif1,h_cd73_hif1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs all GLUT1

popAid = tumor_id; %tumor + ambig
popBid = tumor_id==tumor_id; % all cells
magA = cleandata.CD73;
magB = cleandata.GLUT1;
tic
[spatcorr_cd73_glut1,radfunc_cd73_glut1,parms_cd73_glut1,h_cd73_glut1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs all iNOS

popAid = tumor_id; %tumor + ambig
popBid = tumor_id==tumor_id; % all cells
magA = cleandata.CD73;
magB = cleandata.iNOS;
tic
[spatcorr_cd73_inos,radfunc_cd73_inos,parms_cd73_inos,h_cd73_inos] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs all Arg1

popAid = tumor_id; %tumor + ambig
popBid = tumor_id==tumor_id; % all cells
magA = cleandata.CD73;
magB = cleandata.ARG1;
tic
[spatcorr_cd73_arg1,radfunc_cd73_arg1,parms_cd73_arg1,h_cd73_arg1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% macro cd39-express vs all HIF1-express
popAid = myeloid_id; %macrophages
popBid = myeloid_id==myeloid_id; % all cells
magA = cleandata.CD39;
magB = cleandata.HIF1a; 
tic
[spatcorr_cd39_hif1,radfunc_cd39_hif1,parms_cd39_hif1,h_cd39_hif1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% macro cd39 vs all GLUT1

popAid = myeloid_id; %macrophages
popBid = myeloid_id==myeloid_id; % all cells
magA = cleandata.CD39;
magB = cleandata.GLUT1; 
tic
[spatcorr_cd39_glut1,radfunc_cd39_glut1,parms_cd39_glut1,h_cd39_glut1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% macro cd39 vs all iNOS
popAid = myeloid_id; %macrophages
popBid = myeloid_id==myeloid_id; % all cells
magA = cleandata.CD39;
magB = cleandata.iNOS; 
tic
[spatcorr_cd39_inos,radfunc_cd39_inos,parms_cd39_inos,h_cd39_inos] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% macro cd39 vs all Arg1

popAid = myeloid_id; %macrophages
popBid = myeloid_id==myeloid_id; % all cells
magA = cleandata.CD39;
magB = cleandata.ARG1; 
tic
[spatcorr_cd39_arg1,radfunc_cd39_arg1,parms_cd39_arg1,h_cd39_arg1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs all myeloid
popAid = tumor_id;
popBid = tumor_id==tumor_id;
magA = cleandata.CD73;
magB = myeloid_id;

tic
[spatcorr_cd73_myel,radfunc_cd73_myel,parms_cd73_myel,h_cd73_myel] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs myeloid-CD163
popAid = tumor_id;
popBid = myeloid_id;
magA = cleandata.CD73;
magB = cleandata.CD163;

tic
[spatcorr_cd73_cd163,radfunc_cd73_cd163,parms_cd73_cd163,h_cd73_cd163] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs myeloid-CD11b
popAid = tumor_id;
popBid = myeloid_id;
magA = cleandata.CD73;
magB = cleandata.CD11b;

tic
[spatcorr_cd73_cd11b,radfunc_cd73_cd11b,parms_cd73_cd11b,h_cd73_cd11b] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs myeloid-CD206
popAid = tumor_id;
popBid = myeloid_id;
magA = cleandata.CD73;
magB = cleandata.CD206;

tic
[spatcorr_cd73_cd206,radfunc_cd73_cd206,parms_cd73_cd206,h_cd73_cd206] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc


%% tumor cd73 myeloid-arginase
popAid = tumor_id;
popBid = myeloid_id;
magA = cleandata.CD73;
magB = cleandata.ARG1;

tic
[spatcorr_cd73_myARG1,radfunc_cd73_myARG1,parms_cd73_myARG1,h_cd73_myARG1] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 myeloid-inos
popAid = tumor_id;
popBid = myeloid_id;
magA = cleandata.CD73;
magB = cleandata.iNOS;

tic
[spatcorr_cd73_myNOS,radfunc_cd73_myNOS,parms_cd73_myNOS,h_cd73_myNOS] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73 vs tumor egfr
popAid = tumor_id;
popBid = tumor_id;
magA = cleandata.CD73;
magB = cleandata.EGFR;

tic
[spatcorr_cd73_egfr,radfunc_cd73_egfr,parms_cd73_egfr,h_cd73_egfr] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc
