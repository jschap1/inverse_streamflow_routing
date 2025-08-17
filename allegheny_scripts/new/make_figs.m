% Make figures for the ISR paper
%
% 10/7/2023 JRS

clean
cd('/home/jschap/Dropbox/ISR/inverse_streamflow_routing')
addpath(genpath('./src/'))
load('./allegheny_data/setup/setup-swot-gage.mat');
i1=60; i2=179;
nt=i2-i1+1;
lw=2;
fs=16;
ms=30;

truth = struct();
true_runoff = true_runoff(:,i1:i2);
truth.true_runoff = true_runoff;
truth.total_runoff = true_runoff;
basin.true_runoff = true_runoff;

tv = 1:120;
gi = (k+1):nt-(k+1);

y20_AG = load('./allegheny_data/results/allgage/Y20AG_m0a1L40T5_swot_15.mat');
y20_LM = load('./allegheny_data/results/allgage/Y20LM_m1a1L40T5_swot_15.mat');
ens_AG = load('./allegheny_data/results/allgage/ensAG_m0a1L40T5_swot_15_logup.mat');
ens_LM = load('./allegheny_data/results/allgage/ensLM_m1a1L40T5_swot_15_logup.mat');

% mean prior runoff
ens_AG.meanpriorrunoff = mean(ens_AG.prior_runoff_ens,3);
ens_LM.meanpriorrunoff = mean(ens_LM.prior_runoff_ens,3);

%% Figure 1 - runoff NSE maps

plt_ISR_results_overview(basin, y20_AG.runoff_prior', ...
    y20_AG.post_runoff', truth, tv, gi)
title('Y20 A/G')
plt_ISR_results_overview(basin, ens_AG.meanpriorrunoff, ...
    ens_AG.mean_post_runoff, truth, tv, gi)
title('Ens A/G')
plt_ISR_results_overview(basin, y20_LM.runoff_prior', ...
    y20_LM.post_runoff', truth, tv, gi)
title('Y20 L/M')
plt_ISR_results_overview(basin, ens_LM.meanpriorrunoff, ...
    ens_LM.mean_post_runoff, truth, tv, gi)
title('Ens L/M')

%% Figure 2 - uncertainty ratio boxplots

% calculate uncertainty ratio

top5 = y20_AG.runoff_prior + 1.96*y20_AG.prior_stddev;
bottom5 = y20_AG.runoff_prior - 1.96*y20_AG.prior_stddev;
y20_AG.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = y20_AG.post_runoff + 1.96*y20_AG.post_stddev;
bottom5 = y20_AG.post_runoff - 1.96*y20_AG.post_stddev;
y20_AG.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = y20_LM.runoff_prior + 1.96*y20_LM.prior_stddev;
bottom5 = y20_LM.runoff_prior - 1.96*y20_LM.prior_stddev;
y20_LM.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = y20_LM.post_runoff + 1.96*y20_LM.post_stddev;
bottom5 = y20_LM.post_runoff - 1.96*y20_LM.post_stddev;
y20_LM.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(ens_AG.prior_runoff_ens, 97.5, 3)';
bottom5 = prctile(ens_AG.prior_runoff_ens, 2.5, 3)';
ens_AG.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(ens_AG.post_runoff, 97.5, 3)';
bottom5 = prctile(ens_AG.post_runoff, 2.5, 3)';
ens_AG.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(ens_LM.prior_runoff_ens, 97.5, 3)';
bottom5 = prctile(ens_LM.prior_runoff_ens, 2.5, 3)';
ens_LM.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(ens_LM.post_runoff, 97.5, 3)';
bottom5 = prctile(ens_LM.post_runoff, 2.5, 3)';
ens_LM.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

figure
boxplot([y20_AG.ur95_prior, y20_AG.ur95_post, ens_AG.ur95_post],...
    'labels', {'A/G prior', 'Y20 post', 'EnsISR post'})
title('95% uncertainty ratio')

figure
boxplot([y20_LM.ur95_prior, y20_LM.ur95_post, ens_LM.ur95_post],...
    'labels', {'L/M prior', 'Y20 post', 'EnsISR post'})
title('95% uncertainty ratio')

%% Figure 3 - ER95 boxplots

y20_AG.er95 = myER95_gauss(truth.true_runoff, ...
    y20_AG.post_runoff', y20_AG.post_stddev');
y20_LM.er95 = myER95_gauss(truth.true_runoff, ...
    y20_LM.post_runoff', y20_LM.post_stddev');

ens_AG.er95 = myER(truth.true_runoff, ens_AG.post_runoff, 95);
ens_LM.er95 = myER(truth.true_runoff, ens_LM.post_runoff, 95);

ens_AG.ens_prior_er95 = myER(truth.true_runoff, ens_AG.prior_runoff_ens, 95);
ens_LM.ens_prior_er95 = myER(truth.true_runoff, ens_LM.prior_runoff_ens, 95);
y20_AG.prior_er95 = myER95_gauss(truth.true_runoff, ...
    y20_AG.runoff_prior', y20_AG.prior_stddev');
y20_LM.prior_er95 = myER95_gauss(truth.true_runoff, ...
    y20_LM.runoff_prior', y20_LM.prior_stddev');

figure
boxplot([y20_AG.prior_er95, y20_AG.er95, ens_AG.ens_prior_er95, ens_AG.er95], 'labels', ...
    {'Y20 Prior (A/G)','Y20 Post','Ens Prior (A/G)', 'Ens Post'})
title('95% exceedence ratio')

figure
boxplot([y20_LM.prior_er95, y20_LM.er95, ens_LM.ens_prior_er95, ens_LM.er95], 'labels', ...
    {'Y20 Prior (L/M)','Y20 Post','Ens Prior (L/M)', 'Ens Post'})
title('95% exceedence ratio')

%% Figure 4 - discharge time series

true_discharge = state_model_dumb(true_runoff',HH);
y20_AG.priorQ = state_model_dumb(y20_AG.runoff_prior, HH);
y20_LM.priorQ = state_model_dumb(y20_LM.runoff_prior, HH);
y20_AG.postQ = state_model_dumb(y20_AG.post_runoff, HH);
y20_LM.postQ = state_model_dumb(y20_LM.post_runoff, HH);
ens_AG.priorQ = state_model_dumb(ens_AG.meanpriorrunoff', HH);
ens_LM.priorQ = state_model_dumb(ens_LM.meanpriorrunoff', HH);
ens_AG.postQ = state_model_dumb(ens_AG.mean_post_runoff', HH);
ens_LM.postQ = state_model_dumb(ens_LM.mean_post_runoff', HH);

gagei = 2; % outlet

y20_AG.kge_prior = myKGE(true_discharge(gi,gagei),y20_AG.priorQ(gi,gagei));
y20_LM.kge_prior = myKGE(true_discharge(gi,gagei),y20_LM.priorQ(gi,gagei));
ens_AG.kge_prior = myKGE(true_discharge(gi,gagei),ens_AG.priorQ(gi,gagei));
ens_LM.kge_prior = myKGE(true_discharge(gi,gagei),ens_LM.priorQ(gi,gagei));

% shadedErrorBar(1:365, squeeze(prior_runoff_ens(kk1,:,:))', {@mean, @std}, ...
%     'lineprops', '-r');


subplot(2,2,1) % Y20 AG
plot(y20_AG.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
y20_AG.kge_post = plot_discharge_ts(y20_AG.postQ(gi,gagei), true_discharge(gi,gagei));
title('Y20 A/G')
plot(y20_AG.gage_w_error(gi,gagei), 'k.', 'markersize', ms)
legend('Prior','Truth','Posterior','Observations')
ylim([0,800])
xlim([30,60])
set(gca, 'fontsize', fs)

subplot(2,2,2) % Y20 LM
plot(y20_LM.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
y20_LM.kge_post = plot_discharge_ts(y20_LM.postQ(gi,gagei), true_discharge(gi,gagei));
title('Y20 L/M')
plot(y20_LM.gage_w_error(gi,gagei), 'k.', 'markersize', ms)
ylim([0,800])
xlim([30,60])
set(gca, 'fontsize', fs)
legend('off')

subplot(2,2,3) % ENS AG
plot(ens_AG.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
ens_AG.kge_post = plot_discharge_ts(ens_AG.postQ(gi,gagei), true_discharge(gi,gagei));
title('Ens A/G')
plot(ens_AG.error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', ms)
ylim([0,800])
xlim([30,60])
set(gca, 'fontsize', fs)
legend('off')

subplot(2,2,4) % ENS LM
plot(ens_LM.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
ens_LM.kge_post = plot_discharge_ts(ens_LM.postQ(gi,gagei), true_discharge(gi,gagei));
title('Ens L/M')
plot(ens_LM.error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', ms)
ylim([0,800])
xlim([30,60])
set(gca, 'fontsize', fs)
legend('off')

%% Plot ensemble of runoff

kk=11; % cell number

figure
h1 = plot(y20_AG.runoff_prior(:,kk), 'b*-', 'linewidth', lw);
hold on
plot(y20_AG.runoff_prior(:,kk) + y20_AG.prior_stddev(:,kk), 'b-', 'linewidth', lw)
plot(y20_AG.runoff_prior(:,kk) - y20_AG.prior_stddev(:,kk), 'b-', 'linewidth', lw)
h2 = plot(y20_AG.post_runoff(:,kk), 'r*-', 'linewidth', lw);
plot(y20_AG.post_runoff(:,kk) + y20_AG.post_stddev(:,kk), 'r-', 'linewidth', lw)
plot(y20_AG.post_runoff(:,kk) - y20_AG.post_stddev(:,kk), 'r-', 'linewidth', lw)
h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);
legend([h1,h2,h3], 'Prior','Posterior','Truth')
xlabel('t')
ylabel('runoff (mm/day)')
title('Y20, A/G')
set(gca, 'fontsize', fs)

%% Plot discharge with uncertainty (Y20 AG)

gg=2; % outlet

% calculate discharge uncertainty (error propagation)
% y20_AG.priorQ_std = propagate_runoff_error(y20_AG.prior_stddev, HH, gg);
% y20_LM.priorQ_std = propagate_runoff_error(y20_LM.prior_stddev, HH, gg);
% y20_AG.postQ_std = propagate_runoff_error(y20_AG.post_stddev, HH, gg);
% y20_LM.postQ_std = propagate_runoff_error(y20_LM.post_stddev, HH, gg);

y20_AG.priorQ_std = 0;
y20_LM.priorQ_std = 0;
y20_AG.postQ_std = 0;
y20_LM.postQ_std = 0;

figure
h1=plot(y20_AG.priorQ(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(y20_AG.priorQ(:,gg) + y20_AG.priorQ_std, 'b-', 'linewidth', 1);
plot(y20_AG.priorQ(:,gg) - y20_AG.priorQ_std, 'b-', 'linewidth', 1);
h2=plot(y20_AG.postQ(:,gg), 'r-*', 'linewidth', lw);
plot(y20_AG.postQ(:,gg) + y20_AG.postQ_std, 'r-', 'linewidth', 1);
plot(y20_AG.postQ(:,gg) - y20_AG.postQ_std, 'r-', 'linewidth', 1);
h3=plot(true_discharge(:,gg), 'k', 'linewidth', lw);
h4 = plot(y20_AG.gage_w_error(:,gg), 'k.', 'markersize', ms);
legend([h1,h2,h3,h4],'Prior','Posterior','Truth', 'Obs')
xlabel('t')
xlim([20,60])
ylim([0,800])
title('Y20 AG outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Plot discharge with uncertainty (Y20 LM)

% calculate discharge uncertainty (error propagation)
y20_LM.priorQ_std = propagate_runoff_error(y20_LM.prior_stddev, HH, gg);
y20_LM.postQ_std = propagate_runoff_error(y20_LM.post_stddev, HH, gg);

gg=2; % outlet

figure
h1=plot(y20_LM.priorQ(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(y20_LM.priorQ(:,gg) + y20_LM.priorQ_std, 'b-', 'linewidth', 1);
plot(y20_LM.priorQ(:,gg) - y20_LM.priorQ_std, 'b-', 'linewidth', 1);
h2=plot(y20_LM.postQ(:,gg), 'r-*', 'linewidth', lw);
plot(y20_LM.postQ(:,gg) + y20_LM.postQ_std, 'r-', 'linewidth', 1);
plot(y20_LM.postQ(:,gg) - y20_LM.postQ_std, 'r-', 'linewidth', 1);
h3=plot(true_discharge(:,gg), 'k', 'linewidth', lw);
h4 = plot(y20_LM.gage_w_error(:,gg), 'k.', 'markersize', ms);
legend([h1,h2,h3,h4],'Prior','Posterior','Truth', 'Obs')
xlabel('t')
xlim([20,60])
ylim([0,800])
title('Y20 LM outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Plot discharge with uncertainty (ENS AG)

gg=2; % outlet
ngage = size(HH,1);
M = size(ens_AG.post_runoff,3);

[ens_AG.priorQmean, ens_AG.priorQstd, ~] = ...
    calc_ensemble_discharge_moments(ens_AG.prior_runoff_ens, HH);
[ens_AG.postQmean, ens_AG.postQstd, ~] = ...
    calc_ensemble_discharge_moments(ens_AG.post_runoff, HH);

figure
h1=plot(ens_AG.priorQmean(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(ens_AG.priorQmean(:,gg) + ens_AG.priorQstd(:,gg), 'b-', 'linewidth', 1);
plot(ens_AG.priorQmean(:,gg) - ens_AG.priorQstd(:,gg), 'b-', 'linewidth', 1);
h2=plot(ens_AG.postQmean(:,gg), 'r-*', 'linewidth', lw);
plot(ens_AG.postQmean(:,gg) + ens_AG.postQstd(:,gg), 'r-', 'linewidth', 1);
plot(ens_AG.postQmean(:,gg) - ens_AG.postQstd(:,gg), 'r-', 'linewidth', 1);
h3=plot(true_discharge(:,gg), 'k', 'linewidth', lw);
h4 = plot(ens_AG.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);
legend([h1,h2,h3,h4],'Prior','Posterior','Truth', 'Obs')
xlabel('t')
xlim([20,60])
ylim([0,800])
title('Ens AG outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Plot discharge with uncertainty (ENS LM)

% It might be worth expressing uncertainty in terms of quantiles instead of
% standard deviation..

[ens_LM.priorQmean, ens_LM.priorQstd, ~] = ...
    calc_ensemble_discharge_moments(ens_LM.prior_runoff_ens, HH);
[ens_LM.postQmean, ens_LM.postQstd, ~] = ...
    calc_ensemble_discharge_moments(ens_LM.post_runoff, HH);

figure
h1=plot(ens_LM.priorQmean(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(ens_LM.priorQmean(:,gg) + ens_LM.priorQstd(:,gg), 'b-', 'linewidth', 1);
plot(ens_LM.priorQmean(:,gg) - ens_LM.priorQstd(:,gg), 'b-', 'linewidth', 1);
h2=plot(ens_LM.postQmean(:,gg), 'r-*', 'linewidth', lw);
plot(ens_LM.postQmean(:,gg) + ens_LM.postQstd(:,gg), 'r-', 'linewidth', 1);
plot(ens_LM.postQmean(:,gg) - ens_LM.postQstd(:,gg), 'r-', 'linewidth', 1);
h3=plot(true_discharge(:,gg), 'k', 'linewidth', lw);
h4 = plot(ens_LM.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);
legend([h1,h2,h3,h4],'Prior','Posterior','Truth', 'Obs')
xlabel('t')
xlim([20,60])
ylim([0,800])
title('Ens LM outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)
