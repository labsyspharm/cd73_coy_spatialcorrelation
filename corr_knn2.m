function [corrfunc,rad_approx] = corr_knn2(xydataA,xydataB,magvalA,magvalB,k_val)    
  
    [idx,dist] = knnsearch(xydataA,xydataB,'k',k_val,'NSMethod','kdtree');
    rad_approx = mean(dist); %i.e. the k'th neighbor is on average this distance
    
    magvalA = (magvalA-mean(magvalA))/std(magvalA);
    magvalB = (magvalB-mean(magvalB))/std(magvalB);
    
    magdots = magvalA(idx).*magvalB;
    corrfunc = mean(magdots);
end