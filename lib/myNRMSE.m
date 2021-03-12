function fit = myNRMSE(x, xref)
    
       RMSE = sqrt((sum((xref-x).^2))/length(x));
       NRMSE = RMSE/(max(xref)- min(xref));
       
       fit = 1- NRMSE;
end