allcores = unique(cleandata.core_num_uniq);

%% cell-type IDs are from GMM clustering: 1&3 = Tumor/Glia, 2 = Immune [Lymphoid], 4 = Endothelial, 5 = Immune [Myeloid]

%% tumor cd73-express vs macro cd39-express
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor + ambig
popBid = cleandata.celltype==5; % macrophages
magA = log(cleandata.CD73_sec);
magB = log(cleandata.CD39_mpconj);
tic
[spatcorr_cd73_cd39,radfunc_cd73_cd39,parms_cd73_cd39,h_cd73_cd39] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73-express vs endo presence
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor + ambig
popBid = ~popAid; % nontumor
magA = log(cleandata.CD73_sec);
magB = cleandata.celltype==4; %endothelial
tic
[spatcorr_cd73_endo,radfunc_cd73_endo,parms_cd73_endo,h_cd73_endo] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73-express vs all HIF1-express
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor + ambig
popBid = cleandata.celltype~=0; % all cells
magA = log(cleandata.CD73_sec);
magB = log(cleandata.HIF1a);
tic
[spatcorr_cd73_hif1,radfunc_cd73_hif1,parms_cd73_hif1,h_cd73_hif1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% macro cd39-express vs all HIF1-express
popAid = cleandata.celltype==5; %macrophages
popBid = cleandata.celltype~=0; % all cells
magA = log(cleandata.CD39_mpconj);
magB = log(cleandata.HIF1a); 
tic
[spatcorr_cd39_hif1,radfunc_cd39_hif1,parms_cd39_hif1,h_cd39_hif1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% cd31-presence vs all HIF1-express
popAid = cleandata.celltype~=0; % all cells
popBid = cleandata.celltype~=0; % all cells
magA = cleandata.celltype==4; %endothelial
magB = log(cleandata.HIF1a);
tic
[spatcorr_endo_hif1,radfunc_endo_hif1,parms_endo_hif1,h_endo_hif1] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

%% tumor cd73-express vs cd45 presence and expression of cd3, cd8, cd4, pd1, foxp3, cd24
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor + ambig
popBid = ~popAid; % all non-tumor cells
magA = log(cleandata.CD73_sec);
magB = cleandata.celltype==2; % cd45+
tic
[spatcorr_cd73_cd45,radfunc_cd73_cd45,parms_cd73_cd45,h_cd73_cd45] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc

tic
popBid = cleandata.celltype==2;
imm_mark = {'CD3','CD8A','CD4','PDL1','PD1','FOXP3','CD24'};
for i = 1:numel(imm_mark)
    magB = cleandata{:,imm_mark{i}};
    [spatcorr_cd73_imm{i},radfunc_cd73_imm{i},parms_cd73_imm{i},h_cd73_imm{i}] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
end
toc

%% tumor cd73-express vs macro presence and expression of cd11b, cd14, pd1, pdl1, foxp3
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor + ambig
popBid = ~popAid; % all non-tumor cells
magA = log(cleandata.CD73_sec);
magB = cleandata.celltype==5; % macrophage

[spatcorr_cd73_macro,radfunc_cd73_macro,parms_cd73_macro,h_cd73_macro] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);

tic
popBid = cleandata.celltype==5;
macro_mark = {'CD11b','CD14','PD1','PDL1','FOXP3'};
for i = 1:numel(macro_mark)
    magB = cleandata{:,macro_mark{i}};
    [spatcorr_cd73_macmark{i},radfunc_cd73_macmark{i},parms_cd73_macmark{i},h_cd73_macmark{i},skipflag] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
end
toc