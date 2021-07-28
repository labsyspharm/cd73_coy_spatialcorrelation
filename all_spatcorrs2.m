allcores = unique(cleandata.core_num_uniq);

%% tumor cd73-express vs macro expression of cd68, pu1, cd163, cd11b, cd14
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor + ambig
popBid = cleandata.celltype==5; % macrophages (actually monocytes)
magA = log(cleandata.CD73_sec);

tic
macro_mark2 = {'CD68','PU_1','CD163','CD11b','CD14'};
for i = 1:numel(macro_mark2)
    magB = cleandata{:,macro_mark2{i}};
    [spatcorr_cd73_macmark2{i},radfunc_cd73_macmark2{i},parms_cd73_macmark2{i},h_cd73_macmark2{i},skipflag] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
end
toc

%% macro expression of cd39 vs tumor expression of egfr, p53, ki67, pH2ax, p21, SOX2, OLIG2
popAid = cleandata.celltype==5;
popBid = cleandata.celltype==1 | cleandata.celltype==3; %tumor cells
magA = log(cleandata.CD39_mpconj);

tic
tumormark = {'EGFR','p53','KI67','pH2AX','p21','SOX2','OLIG2'};
for i = 1:numel(tumormark)
    magB = cleandata{:,tumormark{i}};
    [spatcorr_cd39_tumormark{i},radfunc_cd39_tumormark{i},parms_cd39_tumormark{i},h_cd39_tumormark{i},skipflag] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
end
toc


%% tumor expression of cd39 vs all expression of egfr, p53, ki67
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor cells
popBid = cleandata.celltype~=0; %all cells, INCLUDING TUMOR
magA = log(cleandata.CD73_sec);

tic
tumormark = {'EGFR','p53','KI67','pH2AX','p21','SOX2','OLIG2'};
for i = 1:3
    magB = cleandata{:,tumormark{i}};
    [spatcorr_cd73_tumormark{i},radfunc_cd73_tumormark{i},parms_cd73_tumormark{i},h_cd73_tumormark{i},skipflag] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
end
toc

%% all cd73-express vs all cd39-express
popAid = cleandata.celltype~=0; %all cells
popBid = cleandata.celltype~=0; %all cells
magA = log(cleandata.CD73_sec);
magB = log(cleandata.CD39_mpconj);
tic
[spatcorr_all_cd73_cd39,radfunc_all_cd73_cd39,parms_all_cd73_cd39,h_all_cd73_cd39] = ...
    subpop_spatcorr_express_self(cleandata,popAid,popBid,magA,magB,allcores);
toc

