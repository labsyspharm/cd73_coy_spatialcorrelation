function [spatcorrs,radfuncs,parms,h,skipflag] = subpop_spatcorr_express_self(data,popAid,popBid,magA,magB,allcores)
    spatcorrs = cell(numel(allcores),1);
    radfuncs = cell(numel(allcores),1);
    skipflag = zeros(numel(allcores),1);
    h = zeros(numel(allcores),1);
    parms = zeros(numel(allcores),2);

    for i = 1:numel(allcores)
        coreid = data.core_num_uniq == allcores(i);
        if sum(popAid&coreid)<30 || sum(popBid&coreid)<30
            skipflag(i) = 1;
            continue
        end
        xydata = [data.CentroidRow,data.CentroidCol];
        magvalA = magA(coreid&popAid);
        magvalB = magB(coreid&popBid);
        xyA = xydata(coreid&popAid,:);
        xyB = xydata(coreid&popBid,:);
        
        [spatcorrs{i},radfuncs{i}] = corr_knn2(xyA,xyB,magvalA,magvalB,min([sum(popAid&coreid),sum(popBid&coreid),200]));
        
        if(any(isnan(spatcorrs{i}))||any(isnan(radfuncs{i})))
            continue
        end
        
        length = numel(radfuncs{i});
        f = fit(radfuncs{i}(2:length)',spatcorrs{i}(2:length)','exp1','StartPoint',[0,0]);
        fitparms = confint(f,.95);
        [h(i),parm_exp] = cust_hyp(fitparms);
        parms(i,1) = parm_exp(1);
        parms(i,2) = -1/parm_exp(2);
    end
end

function [h,parm_exp] = cust_hyp(parms)
    parm_exp = mean(parms);
    if all(parms(:,1)>0)
        h = 1;
    elseif all(parms(:,1)<0)
        h = -1;
    else
        h = 0;
    end
end