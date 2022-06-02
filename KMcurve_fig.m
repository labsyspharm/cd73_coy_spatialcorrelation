%% Run the three lines below to generate dependent values for the rest of the script
%load('CD73NatComm_Coy2021_CyCIF_TMAdata.mat')
%subtype_comp_fig
%all_spatcorrs

%%
cleanmetadata = metadata(~lostcores,:);
tumorid = cleandata.celltype==1 | cleandata.celltype==3;
macroid = cleandata.celltype==5;
macro_core = grpstats(macroid,cleandata.core_num_uniq)>0;
survival_list = metadata.pfs(~lostcores);
tumorexpress_core = grpstats(log(cleandata.CD73_sec(tumorid)),cleandata.core_num_uniq(tumorid));


cd73cut = 8;
glio_id = ismember(core_diag,glio_diag);
idh1wt_id = cleanmetadata.idh1_cd==2;

figure()
    cd73_cycif = tumorexpress_core<cd73cut; %tumor cd73
    makeKMcurve(survival_list,glio_id,idh1wt_id,cleanmetadata,cd73_cycif)
    title('Tumor CD73')

figure()
    cd73_cycif = h_cd73_cd39~=1; %cd73-cd39 interaction
    makeKMcurve(survival_list,glio_id,idh1wt_id,cleanmetadata,cd73_cycif)
    title('CD73-Tumor and CD39-Macro interaction')

    %saveas(gcf,'KMcurves.pdf')
function makeKMcurve(survival_list,glio_id,idh1wt_id,cleanmetadata,cd73_cycif)
    survival_patient = grpstats(survival_list,cleanmetadata.patient_id);
    status_patient = logical(grpstats(glio_id&cd73_cycif&idh1wt_id,cleanmetadata.patient_id));

    ecdf(survival_patient(status_patient),'function','survivor')
    hold on
    ecdf(survival_patient(~status_patient),'function','survivor')
    legend('CD73-lo','CD73-hi')
    ylabel('PFS')
    xlabel('Days')
    p = ranksum(survival_patient(status_patient),survival_patient(~status_patient));
    text(mean(xlim),mean(ylim),['p=' num2str(p,2)])
end