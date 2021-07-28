load CD73NatComm_Coy2021_CyCIF_TMAdata.mat

allcores = unique(cleandata.core_num_uniq);

%specify two cell subpopulations to compare
popAid = cleandata.celltype==1 | cleandata.celltype==3; %tumor 
popBid = cleandata.celltype==5; % monocytes

%specify which markers to compare for each subpopulation
magA = log(cleandata.CD73_sec);
magB = log(cleandata.CD39_mpconj);

tic
%compute correlation between Tumor CD73 expression and Monocyte CD39 expression 
%as a function of distance for each core

%spatcorr gives the Pearson correlation for the k'th neighbor
%radfunc gives the average inter-cell distance d of the k'th neighbor
%parms gives the exponential fit parameters 'a', 'b' in a*exp(-d/b)
%h indicates if the parameter 'a' has a 95%-confidence interval containing
%0 or not (0 indicates yes, -1 and 1 indicate no)
[spatcorr_cd73_cd39,radfunc_cd73_cd39,parms_cd73_cd39,h_cd73_cd39] = ...
    subpop_spatcorr_express(cleandata,popAid,popBid,magA,magB,allcores);
toc


figure('units','normalized','outerposition',[0 0 1 1]) %example result from one core
subplot(1,2,1) %spatial correlation and fit
    example_coreid = 251; %choice of core to show as example
    
    plot(radfunc_cd73_cd39{example_coreid},spatcorr_cd73_cd39{example_coreid})
    xlim([0 radfunc_cd73_cd39{example_coreid}(end)])
    hold on
    a = parms_cd73_cd39(example_coreid,1); b = parms_cd73_cd39(example_coreid,2);
    a = a*abs(h_cd73_cd39(example_coreid));
    fplot(@(t) a*exp(-t/b),[0 radfunc_cd73_cd39{example_coreid}(end)])
    title({'Tumor-CD73 vs Monocyte-CD39', 'Spatial Correlation'})
    xlabel('Average inter-cell distance (um)'); ylabel('Correlation')
    legend('Data estimate', 'Exponential Fit')
    
    if h_cd73_cd39(example_coreid)==0
        param_text = 'Not Significant (95% CI at 0 um contains 0)';
    else
        param_text = ['Correlation at 0 um is ' num2str(a,2)];
    end
    text(0.5*radfunc_cd73_cd39{example_coreid}(end),0.5*max(spatcorr_cd73_cd39{example_coreid}),param_text,'HorizontalAlignment','center')

subplot(1,2,2) %Spatial coordinates of cells
    coreid = cleandata.core_num_uniq==allcores(example_coreid);
    colorA = zscore(magA)/4+0.5;
    colorB = zscore(magB)/4+0.5;
    colorvec = [popAid.*colorA, popAid~=popAid ,popBid.*colorB] + 0.5*(~popAid & ~popBid)*[1 1 1];
    scatter(cleandata.CentroidRow(coreid),cleandata.CentroidCol(coreid),30,colorvec(coreid,:),'filled','HandleVisibility','off')
    daspect([1 1 1])
    
    xlabel('X-coordinates'); ylabel('Y-coordinates')
    title('Single-cell Coordinates and Expression')
    hold on; scatter(NaN,NaN,1,'r','filled'); 
    scatter(NaN,NaN,'b','filled');
    scatter(NaN,NaN,1,0.5*[1 1 1],'filled')
    legend('Tumor-CD73','Monocyte-CD39','Other cells','Location','southwest')
