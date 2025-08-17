%% How does estimated discharge compare to observations?

prior.mean_runoff = mean(prior_runoff_ens,3);
posterior.mean_runoff = mean(post_runoff,3);

prior.meanQ = state_model_dumb(prior.mean_runoff', HH);
posterior.meanQ = state_model_dumb(posterior.mean_runoff', HH);

gagei = 1;

kge_prior = plot_discharge_ts(prior.meanQ(:,gagei), true_discharge(:,gagei));

figure
plot(prior.meanQ(:,gagei), 'b', 'linewidth', lw)
hold on
kge_post = plot_discharge_ts(posterior.meanQ(:,gagei), true_discharge(:,gagei))
title('Discharge (M=200, 10 days, no err)')
plot(error_corrupted_discharge_meas(:,gagei), 'k.', 'markersize', 20)
% xlim([100,160])
legend('Prior','Truth','Posterior','Observations')
