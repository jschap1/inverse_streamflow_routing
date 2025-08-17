function [kge, rmse, nse]= plt_ISR_results_overview(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, varargin)

if nargin<6
    nt = length(tv);
    gi = 1:nt;
else
    gi = varargin{:};
end

n = size(basin.true_runoff,1);

% figure,plot(mean_posterior_runoff(7,:))
% hold on, plot(truth.total_runoff(7,:))
kge = struct();
rmse = struct();
nse = struct();

figure
subplot(1,3,1)
[nse.prior, kge.prior, rmse.prior, nsemap] = plot_gofmaps(basin, mean_prior_runoff', truth.total_runoff, gi);
hold on
plot(basin.gage_lon, basin.gage_lat, 'r.')
title('Prior mean (NSE)')
% title('Prior mean (KGE)')
subplot(1,3,2)
[nse.posterior, kge.posterior, rmse.posterior, nsemap] = plot_gofmaps(basin, mean_posterior_runoff', truth.total_runoff, gi);
hold on
plot(basin.gage_lon, basin.gage_lat, 'r.')
% title('Posterior mean (KGE)')
title('Posterior mean (NSE)')
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, make_map(basin, nse.posterior-nse.prior), 'sdf')
hold on
plot(basin.gage_lon, basin.gage_lat, 'r.')
title('Improvement in NSE')

% nn = 21; % must be odd
% cmap = cbrewer2('seq','YlGnBu',nn); % blues at bottom
% colormap([1,0,0;cmap]); % black for negative values

figure
plot(tv, mean(mean_prior_runoff,1),'b', 'linewidth', 2)
hold on
plot(tv, mean(mean_posterior_runoff,1),'r', 'linewidth', 2)
plot(tv, mean(truth.total_runoff,1), 'k', 'linewidth', 2)
legend('prior','posterior','truth')
title('Basin avg runoff')
set(gca, 'fontsize', 16)

return
