%load('CD73NatComm_Coy2021_CyCIF_TMAdata.mat')

glio_diag = [1:7];
olig_diag = [8:11];
lostcores = isnan(metadata.core_num_uniq);
core_diag = metadata.diagnosis_cd(~lostcores);
glio_id = ismember(core_diag,glio_diag);
olig_id = ismember(core_diag,olig_diag);
tonsil_id = core_diag==16;

boxplot_id = olig_id+2*glio_id+3*tonsil_id;
nullid_core = ~boxplot_id;
boxnames = {'Olig','Glio','Tonsil (Control)'};

macrocounts = grpstats(cleandata.celltype==5,cleandata.core_num_uniq,'sum'); %total macro
endocounts = grpstats(cleandata.celltype==4,cleandata.core_num_uniq,'sum'); %total endo
immcounts = grpstats(cleandata.celltype==2,cleandata.core_num_uniq,'sum'); %total non-macro immune
corecounts = grpstats(cleandata.QCflag,cleandata.core_num_uniq,'sum'); % total cell counts

figure('units','normalized','outerposition',[0 0 1 0.5])
subplot(1,3,1)
    boxplot(macrocounts(~nullid_core)./corecounts(~nullid_core),boxplot_id(~nullid_core))
    p = ranksum(macrocounts(glio_id)./corecounts(glio_id),macrocounts(olig_id)./corecounts(olig_id));
    set(gca,'xticklabels',boxnames); ylabel('Fraction of cells')
    title('Immune (Monocytic) fraction in cores')
    line([1,2],mean(ylim)*[1,1],'Color','g')
    text(1.5,mean(ylim),['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')
subplot(1,3,2)
    boxplot(endocounts(~nullid_core)./corecounts(~nullid_core),boxplot_id(~nullid_core))
    p = ranksum(endocounts(glio_id)./corecounts(glio_id),endocounts(olig_id)./corecounts(olig_id));
    set(gca,'xticklabels',boxnames); ylabel('Fraction of cells')
    title('Endothelial fraction in cores');
subplot(1,3,3)
    boxplot(immcounts(~nullid_core)./corecounts(~nullid_core),boxplot_id(~nullid_core))
    yyaxis right
    boxplot(immcounts(~nullid_core)./corecounts(~nullid_core),boxplot_id(~nullid_core))
    yyaxis left; ylim([0 0.2]); ylabel('Fraction of cells')
    yyaxis right; ylim([0.1 1]); ylabel('Fraction of cells')
    line(2.5*[1 1],ylim)
    p = ranksum(immcounts(glio_id)./corecounts(glio_id),immcounts(olig_id)./corecounts(olig_id));
    set(gca,'xticklabels',boxnames);
    title('Immune (Non-monocytic) fraction in cores');
    line([1,2],mean(ylim)*[1,1],'Color','g')
    text(1.5,mean(ylim),['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')

    
saveas(gcf,'Core_compositions.png')