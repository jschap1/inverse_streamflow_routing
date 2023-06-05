function plt_ISR_results_overview(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, varargin)

if nargin<6
    nt = length(tv);
    gi = 1:nt;
else
    gi = varargin{:};
end

n = size(basin.true_runoff,1);

figure
subplot(1,2,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, mean_prior_runoff', truth.total_runoff, gi);
title('Prior mean (NSE)')
% title('Prior mean (KGE)')
subplot(1,2,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, mean_posterior_runoff', truth.total_runoff, gi);
% title('Posterior mean (KGE)')
title('Posterior mean (NSE)')

figure
plot(tv, mean(mean_prior_runoff,1),'b', 'linewidth', 2)
hold on
plot(tv, mean(mean_posterior_runoff,1),'r', 'linewidth', 2)
plot(tv, mean(truth.total_runoff,1), 'k', 'linewidth', 2)
legend('prior','posterior','truth')
title('Basin avg runoff')
set(gca, 'fontsize', 16)

return
