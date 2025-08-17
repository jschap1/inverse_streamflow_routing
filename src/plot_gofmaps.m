% Plot GOF maps
%
% Plots maps of goodness-of-fit metrics for ISR
%
% INPUTS
% isr = ISR runoff estimate (n by nt)
% truth = true runoff (n by nt)

function [nse, kge, rmse, nsemap] = plot_gofmaps(basin, isr, truth, gi)


isr = isr(gi,:);
truth = truth(:,gi);

[nt,n] = size(isr);

rmse = zeros(n,1);
kge = zeros(n,1);
nse = zeros(n,1);
for kk=1:n
%     rmse(kk) = myRMSE(truth.totR(reorder(kk), :)', isr.r(:,kk));
%     nse(kk) = myNSE(truth.totR(reorder(kk), :)', isr.r(:,kk));
%     kge(kk) = myKGE(truth.totR(reorder(kk), :)', isr.r(:,kk));
    
    true_runoff = truth(kk,:)';
    
    if all(isnan(true_runoff))
        rmse(kk) = NaN;
        nse(kk) = NaN;
        kge(kk) = NaN;
    else
        rmse(kk) = myRMSE(true_runoff, isr(:,kk));
        nse(kk) = myNSE(true_runoff, isr(:,kk));
        kge(kk) = myKGE(true_runoff, isr(:,kk));
    end
    
%     figure
%     plot(true_runoff)
%     hold on
%     plot(isr(:,kk))
%     legend('truth','estimate')
    
end

nsemap = make_map(basin, nse);
kgemap = make_map(basin, kge);
rmsemap = make_map(basin, rmse);

%figure
% plotraster(basin.lonv, basin.latv, kgemap, 'KGE');
% colormap(bluewhitered(256)), colorbar
% caxis([0,1])

% figure
plotraster(basin.lonv, basin.latv, nsemap, 'NSE');
% colormap(bluewhitered(256)), colorbar
caxis([0,1])

% figure
% plotraster(basin.lonv, basin.latv, rmsemap, 'RMSE');
% % colormap(bluewhitered(256)), colorbar

% figure
% plot(rmse)
% hold on
% plot(kge)
% plot(nse)
% legend('rmse','kge','nse')

% Min RMSE does not map to max NSE or KGE, oddly

% [~, maxi] = max(kge)

% disp(['mean NSE: ' num2str(nanmean(nse))])
disp(['median NSE: ' num2str(nanmedian(nse))])

% disp(['min KGE: ' num2str(nanmin(kge))])
% disp(['max KGE: ' num2str(nanmax(kge))])
% disp(['mean KGE: ' num2str(nanmean(kge))])
% disp(['median KGE: ' num2str(nanmedian(kge))])
% 
% disp(['mean RMSE: ' num2str(nanmean(rmse))])
% disp(['median RMSE: ' num2str(nanmedian(rmse))])

return