%plot of Tumor-CD73 relations
figure(); title('CHOP TMA Tumor-CD73 Spatial Correlations (N=183 Cores)')
hold on
plotdensity(parms_cd73_cd39(~ctrl_core_id,1).*abs(h_cd73_cd39(~ctrl_core_id)),'r','Myeloid CD39')
plotdensity(parms_cd73_cd11b(~ctrl_core_id,1).*abs(h_cd73_cd11b(~ctrl_core_id)),'b','Myeloid CD11b')
plotdensity(parms_cd73_cd163(~ctrl_core_id,1).*abs(h_cd73_cd163(~ctrl_core_id)),'m','Myeloid CD163')
plotdensity(parms_cd73_cd206(~ctrl_core_id,1).*abs(h_cd73_cd206(~ctrl_core_id)),'g','Myeloid CD206')
plotdensity(parms_cd73_myARG1(~ctrl_core_id,1).*abs(h_cd73_myARG1(~ctrl_core_id)),'c','Myeloid ARG1')
plotdensity(parms_cd73_myNOS(~ctrl_core_id,1).*abs(h_cd73_myNOS(~ctrl_core_id)),'k','Myeloid iNOS')
plotdensity(parms_cd73_egfr(~ctrl_core_id,1).*abs(h_cd73_egfr(~ctrl_core_id)),'y','Tumor EGFR')
xlim([-0.5 1]); xlabel('Correlation Strength')
ylabel('Proportion of Cores')
legend()

%plot of HIF1 relations
figure(); title('CHOP TMA Tumor-CD73 Hypoxia Spatial Correlations (N=183 Cores)')
hold on
plotdensity(parms_cd73_hif1(~ctrl_core_id,1).*abs(h_cd73_hif1(~ctrl_core_id)),'r','All-HIF1a')
plotdensity(parms_cd73_glut1(~ctrl_core_id,1).*abs(h_cd73_glut1(~ctrl_core_id)),'g','All-Glut1')
plotdensity(parms_cd73_inos(~ctrl_core_id,1).*abs(h_cd73_inos(~ctrl_core_id)),'b','All-iNOS')
plotdensity(parms_cd73_arg1(~ctrl_core_id,1).*abs(h_cd73_arg1(~ctrl_core_id)),'k','All-ARG1')
xlim([-0.5 1]); xlabel('Correlation Strength')
ylabel('Proportion of Cores')
legend()

figure(); title('CHOP TMA Myeloid-CD39 Hypoxia Spatial Correlations (N=183 Cores)')
hold on
plotdensity(parms_cd39_hif1(~ctrl_core_id,1).*abs(h_cd39_hif1(~ctrl_core_id)),'r','All-HIF1a')
plotdensity(parms_cd39_glut1(~ctrl_core_id,1).*abs(h_cd39_glut1(~ctrl_core_id)),'g','All-Glut1')
plotdensity(parms_cd39_inos(~ctrl_core_id,1).*abs(h_cd39_inos(~ctrl_core_id)),'b','All-iNOS')
plotdensity(parms_cd39_arg1(~ctrl_core_id,1).*abs(h_cd39_arg1(~ctrl_core_id)),'k','All-ARG1')
xlim([-0.5 1]); xlabel('Correlation Strength')
ylabel('Proportion of Cores')
legend()


function plotdensity(data,color,lgdlabel)
    bw = 0.05;
    [f,x] = ksdensity(data(~isnan(data)),linspace(-0.5,1,300),'Bandwidth',bw);
    [~,pval] = ttest(data(~isnan(data)));
    plot(x,f,color,'LineWidth',1,'DisplayName',[lgdlabel ' p=' num2str(pval,'%.2e') ' rho~' num2str(mean(data),'%.2f')])
end