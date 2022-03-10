load('CD73NatComm_Coy2021_CHOPTMAdata.mat')
table_cd73cd39 = grpstats(cleandata{:,{'Core_Map_ID','CD73','CD39'}},cleandata.Core_Map_ID);
table_cd73cd39 = array2table(table_cd73cd39,'VariableNames',{'Core_Map_ID','CD73','CD39'});

table_cd73cd39(table_cd73cd39.Core_Map_ID==0,:) = [];

for i = 1:size(table_cd73cd39,1)
    j = find(mutationdata.Core_Map_ID==table_cd73cd39.Core_Map_ID(i));
    table_cd73cd39.CD73index(i) = mutationdata.CD73Index(j);
    table_cd73cd39.CD73prot(i) = mutationdata.NT5EProtZ(j);
    table_cd73cd39.CD73rna(i) = mutationdata.NT5ERNAZ(j);
    table_cd73cd39.CD39prot(i) = mutationdata.ENTPD1ProtZ(j);
    table_cd73cd39.CD39rna(i) = mutationdata.ENTPD1RNAZ(j);
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
    scatter(table_cd73cd39.CD73,table_cd73cd39.CD73prot,'filled')
    [rho,pval] = corr(table_cd73cd39.CD73,table_cd73cd39.CD73prot,'Rows','complete');
    xlabel('CyCIF mean patient CD73'); ylabel('CD39 Protein Z-score')
    title(['Pearson Corr = ' num2str(rho) '; p-val = ' num2str(pval)])
    axis square
subplot(2,2,2)
    scatter(table_cd73cd39.CD73,table_cd73cd39.CD73rna,'filled')
    [rho,pval] = corr(table_cd73cd39.CD73,table_cd73cd39.CD73rna,'Rows','complete');
    xlabel('CyCIF mean patient CD73'); ylabel('CD39 RNA Z-score')
    title(['Pearson Corr = ' num2str(rho) '; p-val = ' num2str(pval)])
    axis square
subplot(2,2,3)
    scatter(table_cd73cd39.CD39,table_cd73cd39.CD39prot,'filled')
    [rho,pval] = corr(table_cd73cd39.CD39,table_cd73cd39.CD39prot,'Rows','complete');
    xlabel('CyCIF mean patient CD39'); ylabel('CD39 Protein Z-score')
    title(['Pearson Corr = ' num2str(rho) '; p-val = ' num2str(pval)])
    axis square
subplot(2,2,4)
    scatter(table_cd73cd39.CD39,table_cd73cd39.CD39rna,'filled')
    [rho,pval] = corr(table_cd73cd39.CD39,table_cd73cd39.CD39rna,'Rows','complete');
    xlabel('CyCIF mean patient CD39'); ylabel('CD39 RNA Z-score')
    title(['Pearson Corr = ' num2str(rho) '; p-val = ' num2str(pval)])
    axis square
saveas(gcf,'cd73cd39_choptma_rnaprotcycif.pdf')
saveas(gcf,'cd73cd39_choptma_rnaprotcycif.png')
figure()
    boxplot(table_cd73cd39.CD73,table_cd73cd39.CD73index)
    xlabel('CD73 Index'); ylabel('CyCIF CD73 Mean Expression')
    [rho,pval] = corr(table_cd73cd39.CD73,table_cd73cd39.CD73index,'Rows','Complete');
    title(['Pearson Corr = ' num2str(rho) '; p-val = ' num2str(pval)])
saveas(gcf,'cd73_choptma_indexcycif.pdf')