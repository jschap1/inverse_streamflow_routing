% Make figures for the ISR paper
%
% 11/6/2023 JRS
% This version is for the Y20 ensemble case

clear
cd('/home/jschap/Dropbox/ISR/inverse_streamflow_routing/allegheny_data/results/')
addpath(genpath('./src/'))
load('../setup/setup-swot-gage.mat');
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

% Daily measurements at outlet, 5% meas error (Cases 1-4)
y20_AG = load('./outlet_only/Cases_1-4/Case2_Y20LM_m1a1L40T5_daily_5_unbiased_M1-500.mat'); % Case 2 LM errors, Y20 prop
% y20_AG = load('./outlet_only/Cases_1-4/Case1_Y20AG_m0a1L40T5_daily_5_unbiased_M1-500.mat'); % Case 1 AG errors, Y20 prop
% y20_LM = load('./outlet_only/Cases_1-4/Case2_Y20LM_m1a1L40T5_daily_5_unbias_1-500.mat'); % Case 2 LM errors, Y20 prop
y20_LM = load('./outlet_only/Cases_1-4/Case3.1_Y20LM_m0a1L40T5_daily_5_unbiased_flat_alpha4_M1-500.mat'); % Case 3.1 LM errors, Y20 flat
% y20_LM = load('./outlet_only/Cases_1-4/Case3.2_Y20LM_m0a1L40T5_daily_5_unbiased_flat_alpha0.25_M1-500.mat'); % Case 3.2 LM errors, Y20 flat
% y20_TMPA = load('./outlet_only/Cases_1-4/Case4_Y20AG_m0a1L40T5_daily_5_tmpa_outlet_only.mat')

% AG vs. LM, Y20 vs. Ens
% Unbiased prior, 15% error, 10-day meas, 14 gauges
ens_LM = load('./allgage/Ens/ensLM_m1a1L40T5_swot_15.mat');
ens_AG = ens_LM;
% ens_AG = load('./allegheny_data/results/outlet_only/ensLM_m1a1L40T5_swot_nomeaserr.mat');
% y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_unbias_1-500.mat');
% y20_LM = load('./allegheny_data/results/allgage/Y20LM_m1a1L40T5_swot_15_unbias_flat_1-M.mat');

% Y20LM_m1a1L40T5_daily_5_tmpa.mat
% y20_LM = load('./allegheny_data/results/outlet_only/Y20LM_m1a1L40T5_daily_5_unbias_1-500.mat');
% y20_LM = load('./allegheny_data/results/allgage/Y20LM_m1a1L40T5_daily_5_unbias_flat_1-500');

% flat AG case
% y20_AG = load('./allegheny_data/results/allgage/Y20AGflat_m0L40T5_swot_15_1-M.mat');
% tmp = load('allegheny_data/errors/Y20AGflat_m0L40T5_prior_ensemble.mat');
% y20_AG.prior_runoff = permute(tmp.prior_runoff_ens, [2, 1, 3]); % needed for AGflat runoff prior...
%ens_AG = load('./allegheny_data/results/allgage/ensAGflat_m0a1L40T5_swot_15.mat');

% y20_LM = load('./allegheny_data/results/allgage/Y20LM_m1a1L40T5_swot_15_1-M.mat');
% ens_AG = load('./allegheny_data/results/allgage/ensAG_m0a1L40T5_swot_15.mat');
% ens_LM = load('./allegheny_data/results/allgage/ensLM_m1a1L40T5_swot_15_logup.mat');

%%

y20_AG.prior_runoff = permute(y20_AG.prior_runoff_ens, [1,2,3]); % if nec
y20_LM.prior_runoff = permute(y20_LM.prior_runoff_ens, [1,2,3]); % if nec

% mean prior runoff
ens_AG.meanpriorrunoff = mean(ens_AG.prior_runoff_ens,3);
ens_LM.meanpriorrunoff = mean(ens_LM.prior_runoff_ens,3);

y20_AG.mean_post_runoff = mean(y20_AG.post_runoff, 3);
y20_AG.mean_prior_runoff = mean(y20_AG.prior_runoff, 3); 

y20_LM.mean_post_runoff = mean(y20_LM.post_runoff, 3);
y20_LM.mean_prior_runoff = mean(y20_LM.prior_runoff, 3)'; 

% Set up colors for figures
colors = [0,0,255; % dark blue (Ens prior)
           0, 255, 0; % green (Y20 posterior)
           255,0,0;  % red (Ens posterior)
           0,255,255]/255; % cyan (Y20 prior)
       
%% Reduction in runoff uncertainty (overview)

mprs = mean(y20_AG.prior_stddev,3);
mpts = mean(y20_AG.post_stddev,3);
figure,imagesc(mprs - mpts), colorbar, title('Y20 AG')

mprs = mean(y20_LM.prior_stddev,3);
mpts = mean(y20_LM.post_stddev,3);
figure,imagesc(mprs - mpts), colorbar, title('Y20 LM')
       
%% Figure 1 - runoff NSE maps

plt_ISR_results_overview(basin, y20_AG.mean_prior_runoff, ...
    y20_AG.mean_post_runoff', truth, tv, gi)
title('Y20 A/G')

plt_ISR_results_overview(basin, ens_AG.meanpriorrunoff, ...
    ens_AG.mean_post_runoff, truth, tv, gi)
title('Ens A/G')

plt_ISR_results_overview(basin, y20_LM.mean_prior_runoff, ...
    y20_LM.mean_post_runoff', truth, tv, gi)
title('Y20 L/M')

plt_ISR_results_overview(basin, ens_LM.meanpriorrunoff, ...
    ens_LM.mean_post_runoff, truth, tv, gi)
title('Ens L/M')

% plt_ISR_results_overview(basin, ens_AG.prior_runoff_std', ...
%     ens_AG.post_runoff_std', truth, tv, gi)
% title('Y20 L/M')

%% Figure 2 - uncertainty ratio boxplots

% calculate uncertainty ratio

top5 = prctile(y20_AG.prior_runoff, 97.5, 3);
bottom5 = prctile(y20_AG.prior_runoff, 2.5, 3);
y20_AG.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(y20_AG.post_runoff, 97.5, 3);
bottom5 = prctile(y20_AG.post_runoff, 2.5, 3);
y20_AG.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(y20_LM.prior_runoff, 97.5, 3);
bottom5 = prctile(y20_LM.prior_runoff, 2.5, 3);
y20_LM.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

top5 = prctile(y20_LM.post_runoff, 97.5, 3);
bottom5 = prctile(y20_LM.post_runoff, 2.5, 3);
y20_LM.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

% top5 = y20_AG.runoff_prior + 1.96*y20_AG.prior_stddev;
% bottom5 = y20_AG.runoff_prior - 1.96*y20_AG.prior_stddev;
% y20_AG.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);
% 
% top5 = y20_AG.post_runoff + 1.96*y20_AG.post_stddev;
% bottom5 = y20_AG.post_runoff - 1.96*y20_AG.post_stddev;
% y20_AG.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);
% 
% top5 = y20_LM.runoff_prior + 1.96*y20_LM.prior_stddev;
% bottom5 = y20_LM.runoff_prior - 1.96*y20_LM.prior_stddev;
% y20_LM.ur95_prior = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);
% 
% top5 = y20_LM.post_runoff + 1.96*y20_LM.post_stddev;
% bottom5 = y20_LM.post_runoff - 1.96*y20_LM.post_stddev;
% y20_LM.ur95_post = sum(top5 - bottom5,1)'./sum(truth.true_runoff,2);

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
    'labels', {'A/G prior', 'Y20 post', 'EnsISR post'}, 'notch', 'on')
title('AG 95% uncertainty ratio')

figure
boxplot([y20_LM.ur95_prior, y20_LM.ur95_post, ens_LM.ur95_post],...
    'labels', {'L/M prior', 'Y20 post', 'EnsISR post'}, 'notch', 'on')
title('LM 95% uncertainty ratio')

%% Figure 3 - ER95 boxplots

% y20_AG.er95 = myER95_gauss(truth.true_runoff, ...
%     mean(y20_AG.post_runoff,3)', mean(y20_AG.post_stddev,3)');
% y20_LM.er95 = myER95_gauss(truth.true_runoff, ...
%     y20_LM.post_runoff', y20_LM.post_stddev');

y20_AG.er95 = myER(truth.true_runoff, permute(y20_AG.post_runoff, [2,1,3]), 95);
y20_LM.er95 = myER(truth.true_runoff, permute(y20_LM.post_runoff, [2,1,3]), 95);

ens_AG.er95 = myER(truth.true_runoff, ens_AG.post_runoff, 95);
ens_LM.er95 = myER(truth.true_runoff, ens_LM.post_runoff, 95);

y20_AG.prior_er95 = myER(truth.true_runoff, permute(y20_AG.prior_runoff, [2,1,3]), 95);
y20_LM.prior_er95 = myER(truth.true_runoff, permute(y20_LM.prior_runoff, [2,1,3]), 95);

ens_AG.ens_prior_er95 = myER(truth.true_runoff, ens_AG.prior_runoff_ens, 95);
ens_LM.ens_prior_er95 = myER(truth.true_runoff, ens_LM.prior_runoff_ens, 95);

% y20_AG.prior_er95 = myER95_gauss(truth.true_runoff, ...
%     y20_AG.runoff_prior', y20_AG.prior_stddev');
% y20_LM.prior_er95 = myER95_gauss(truth.true_runoff, ...
%     y20_LM.runoff_prior', y20_LM.prior_stddev');

figure
boxplot([y20_AG.prior_er95, y20_AG.er95, ens_AG.ens_prior_er95, ens_AG.er95], 'labels', ...
    {'Y20 Prior (A/G)','Y20 Post','Ens Prior (A/G)', 'Ens Post'})
title('AG 95% exceedence ratio')

figure
boxplot([y20_LM.prior_er95, y20_LM.er95, ens_LM.ens_prior_er95, ens_LM.er95], 'labels', ...
    {'Y20 Prior (L/M)','Y20 Post','Ens Prior (L/M)', 'Ens Post'})
title('LM 95% exceedence ratio')

%% Was the consistent case better than the inconsistent case?

% For LM errors, was Y20 with proportional or flat error assumption better?
% Comparing Case 2 and Case 3. Using y20_AG as name for Case 2 here.

% true_runoff = true_runoff';
% figure,histogram(y20_AG.post_runoff - true_runoff)
% figure,histogram(y20_LM.post_runoff - true_runoff)

nbins = 20; 

RMSE2 = zeros(M,1);
for mm=1:M
    aa = true_runoff(gi,:);
    bb = y20_AG.post_runoff(gi,:,mm);
    RMSE2(mm) = myRMSE(aa(:), bb(:));
end
figure
subplot(2,1,1)
histogram(RMSE2,nbins), title('Case 2 RMSE (LM errors with proportional AG error model)')
hold on
plot([mean(RMSE2) mean(RMSE2)],[0 200], 'k')
text(2,60,['mean RMSE = ', num2str(mean(RMSE2))])
xlim([1,4])

RMSE3 = zeros(M,1);
for mm=1:M
    aa = true_runoff(gi,:);
    bb = y20_LM.post_runoff(gi,:,mm);
    RMSE3(mm) = myRMSE(aa(:), bb(:));
end
subplot(2,1,2)
histogram(RMSE3,nbins)
title('Case 3 RMSE (LM errors with flat AG error model)')
hold on
plot([mean(RMSE3) mean(RMSE3)],[0 200], 'k')
text(2,60,['mean RMSE = ', num2str(mean(RMSE3))])
xlim([1,4])

% compare the two samples to see if RMSE is significantly different
h = ttest2(RMSE2, RMSE3)

figure
boxplot([RMSE2, RMSE3], 'notch', 'on') 

%% Figure 4 - discharge time series

true_discharge = state_model_dumb(true_runoff,HH);
[nt,n,M] = size(y20_LM.post_runoff);

% Calculate discharge for each realization of Y20
for mm=1:M
    y20_AG.priorQ(:,:,mm) = state_model_dumb(y20_AG.prior_runoff(:,:,mm), HH);
    y20_LM.priorQ(:,:,mm) = state_model_dumb(y20_LM.prior_runoff(:,:,mm), HH);
    y20_AG.postQ(:,:,mm) = state_model_dumb(y20_AG.post_runoff(:,:,mm), HH);
    y20_LM.postQ(:,:,mm) = state_model_dumb(y20_LM.post_runoff(:,:,mm), HH);
end

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
plot(mean(y20_AG.priorQ(gi,gagei,:),3), 'b-*', 'linewidth', lw)
hold on
y20_AG.kge_post = plot_discharge_ts(mean(y20_AG.postQ(gi,gagei,:),3), true_discharge(gi,gagei));
title('Y20 A/G')
plot(y20_AG.gage_w_error(gi,gagei), 'k.', 'markersize', ms)
legend('Prior','Truth','Posterior','Observations')
ylim([0,800])
xlim([30,60])
set(gca, 'fontsize', fs)

subplot(2,2,2) % Y20 LM
plot(mean(y20_LM.priorQ(gi,gagei,:),3), 'b-*', 'linewidth', lw)
hold on
y20_LM.kge_post = plot_discharge_ts(mean(y20_LM.postQ(gi,gagei,:),3), true_discharge(gi,gagei));
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

y20_AG.mean_prior_std = mean(y20_AG.prior_stddev, 3);
y20_AG.mean_post_std = mean(y20_AG.post_stddev, 3);
y20_LM.mean_prior_std = mean(y20_LM.prior_stddev, 3);
y20_LM.mean_post_std = mean(y20_LM.post_stddev, 3);

kk=40; % cell number

figure
h1 = plot(y20_AG.mean_prior_runoff(:,kk), 'b*-', 'linewidth', lw);
hold on
plot(y20_AG.mean_prior_runoff(:,kk) + y20_AG.mean_prior_std(:,kk), 'b-', 'linewidth', lw)
plot(y20_AG.mean_prior_runoff(:,kk) - y20_AG.mean_prior_std(:,kk), 'b-', 'linewidth', lw)
h2 = plot(y20_AG.mean_post_runoff(:,kk), 'r*-', 'linewidth', lw);
plot(y20_AG.mean_post_runoff(:,kk) + y20_AG.mean_post_std(:,kk), 'r-', 'linewidth', lw)
plot(y20_AG.mean_post_runoff(:,kk) - y20_AG.mean_post_std(:,kk), 'r-', 'linewidth', lw)
h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);
legend([h1,h2,h3], 'Prior','Posterior','Truth')
xlabel('t')
ylabel('runoff (mm/day)')
title('Y20, A/G')
set(gca, 'fontsize', fs)

%% Plot discharge with uncertainty (Y20 AG)

% should show IQR, rather than stddev... (but first use stddev for
% consistency for report to Steve)

gg=2; % outlet

% calculate discharge uncertainty (error propagation)
% y20_AG.priorQ_std = propagate_runoff_error(y20_AG.prior_stddev, HH, gg);
% y20_LM.priorQ_std = propagate_runoff_error(y20_LM.prior_stddev, HH, gg);
% y20_AG.postQ_std = propagate_runoff_error(y20_AG.post_stddev, HH, gg);
% y20_LM.postQ_std = propagate_runoff_error(y20_LM.post_stddev, HH, gg);

y20_AG.priorQ_mean = mean(y20_AG.priorQ,3);
y20_AG.postQ_mean = mean(y20_AG.postQ,3);
y20_AG.priorQ_std = std(y20_AG.priorQ,[],3);
y20_AG.postQ_std = std(y20_AG.postQ,[],3);

figure
h1=plot(y20_AG.priorQ_mean(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(y20_AG.priorQ_mean(:,gg) + y20_AG.priorQ_std(:,gg), 'b-', 'linewidth', 1);
plot(y20_AG.priorQ_mean(:,gg) - y20_AG.priorQ_std(:,gg), 'b-', 'linewidth', 1);
h2=plot(y20_AG.postQ_mean(:,gg), 'r-*', 'linewidth', lw);
plot(y20_AG.postQ_mean(:,gg) + y20_AG.postQ_std(:,gg), 'r-', 'linewidth', 1);
plot(y20_AG.postQ_mean(:,gg) - y20_AG.postQ_std(:,gg), 'r-', 'linewidth', 1);
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
%y20_LM.priorQ_std = propagate_runoff_error(y20_LM.prior_stddev, HH, gg);
%y20_LM.postQ_std = propagate_runoff_error(y20_LM.post_stddev, HH, gg);

y20_LM.priorQ_mean = mean(y20_LM.priorQ,3);
y20_LM.postQ_mean = mean(y20_LM.postQ,3);
y20_LM.priorQ_std = std(y20_LM.priorQ,[],3);
y20_LM.postQ_std = std(y20_LM.postQ,[],3);

gg=2; % outlet

figure
h1=plot(y20_LM.priorQ_mean(:,gg), 'b-*', 'linewidth', lw);
hold on
plot(y20_LM.priorQ_mean(:,gg) + y20_LM.priorQ_std(:,gg), 'b-', 'linewidth', 1);
plot(y20_LM.priorQ_mean(:,gg) - y20_LM.priorQ_std(:,gg), 'b-', 'linewidth', 1);
h2=plot(y20_LM.postQ_mean(:,gg), 'r-*', 'linewidth', lw);
plot(y20_LM.postQ_mean(:,gg) + y20_LM.postQ_std(:,gg), 'r-', 'linewidth', 1);
plot(y20_LM.postQ_mean(:,gg) - y20_LM.postQ_std(:,gg), 'r-', 'linewidth', 1);
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

%% Are there subtle differences we can tease out among the four cases?

% In particular, we might expect Y20 to perform worse than EnsISR for the
% lognormal error case. Discharge seems largely similar, but perhaps runoff
% is different. Also, is there a difference in terms of how many negative
% runoff values there are?

%% Plot y20_AG runoff estimates with shaded error bar

kk=4;

figure(1)

subplot(1,2,1)

meas_times = find(~isnan(y20_AG.gage_w_error(gi,2))); 

% Show median/mean plus or minus IQR

errbars_prior = 1.96*mean(y20_AG.prior_stddev(gi,kk,:),3); % 
errbars_post = 1.96*mean(y20_AG.post_stddev(gi,kk,:),3); % 

h1 = shadedErrorBar(tv(gi),mean(y20_AG.prior_runoff(gi,kk,:),3),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),mean(y20_AG.post_runoff(gi,kk,:),3),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',mean(y20_AG.prior_runoff(gi,kk,:),3))
nse2 = myNSE(true_runoff(kk,gi)',mean(y20_AG.post_runoff(gi,kk,:),3))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

% plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Runoff at cell ' num2str(kk) ' (Y20 AG Prop)'])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

ylim([-6,20])
grid on
set(gca, 'fontsize', fs)
% shows 95% confidence bounds

%% Plot y20_LM runoff estimates with shaded error bar

figure(1)
subplot(1,2,2)

% Show median/mean plus or minus IQR

errbars_prior = 1.96*mean(y20_LM.prior_stddev(gi,kk,:),3); % 
errbars_post = 1.96*mean(y20_LM.post_stddev(gi,kk,:),3); % 

h1 = shadedErrorBar(tv(gi),mean(y20_LM.prior_runoff(gi,kk,:),3),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),mean(y20_LM.post_runoff(gi,kk,:),3),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',mean(y20_LM.prior_runoff(gi,kk,:),3))
nse2 = myNSE(true_runoff(kk,gi)',mean(y20_LM.post_runoff(gi,kk,:),3))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

% plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Runoff at cell ' num2str(kk) ' (Y20 LM)'])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

ylim([-6,20])
grid on
set(gca, 'fontsize', fs)
% shows 95% confidence bounds

%% Plot ensemble of runoff (Ens_AG) with shadederrorbar

% Show median +- IQR
ens_AG.prior_runoff_med = median(ens_AG.prior_runoff_ens,3)';
ens_AG.post_runoff_med = median(ens_AG.post_runoff,3)';
ens_AG.prior_runoff_IQR = iqr(ens_AG.prior_runoff_ens,3)';
ens_AG.post_runoff_IQR = iqr(ens_AG.post_runoff,3)';
ens_AG.prior_runoff_std = std(ens_AG.prior_runoff_ens, [], 3)';
ens_AG.post_runoff_std = std(ens_AG.post_runoff, [], 3)';
% Better yet, show median, 2.5th percentile, and 97.5th percentile
% to reflect UR95 better and avoid making it look like there can be negative values 

ens_AG.prior975 = prctile(ens_AG.prior_runoff_ens, 97.5, 3)';
ens_AG.prior25 = prctile(ens_AG.prior_runoff_ens, 2.5, 3)';
ens_AG.post975 = prctile(ens_AG.post_runoff, 97.5, 3)';
ens_AG.post25 = prctile(ens_AG.post_runoff, 2.5, 3)';

% errbars_prior = [ens_AG.prior975(:,kk), ens_AG.prior25(:,kk)];
% errbars_post = [ens_AG.post975(:,kk), ens_AG.post25(:,kk)];

% errbars_prior = ens_AG.prior_runoff_IQR(:,kk);
% errbars_post = ens_AG.post_runoff_IQR(:,kk);

ens_AG.meanpriorrunoff = mean(ens_AG.prior_runoff_ens,3);
ens_LM.meanpriorrunoff = mean(ens_LM.prior_runoff_ens,3);

figure(1)
subplot(2,2,3)
kk=40;

errbars_prior = 1.96*ens_AG.prior_runoff_std(gi,kk);
errbars_post = 1.96*ens_AG.post_runoff_std(gi,kk);

h1 = shadedErrorBar(tv(gi),ens_AG.meanpriorrunoff(kk,gi),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),ens_AG.mean_post_runoff(kk,gi),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

title(['Runoff at cell ' num2str(kk) ' (Ens AG)'])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
ylim([-6,20])
grid on
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,gi)', ens_AG.meanpriorrunoff(kk,gi)')
myKGE(true_runoff(kk,gi)', ens_AG.mean_post_runoff(kk,gi)')

myNSE(true_runoff(kk,gi)', ens_AG.meanpriorrunoff(kk,gi)')
myNSE(true_runoff(kk,gi)', ens_AG.mean_post_runoff(kk,gi)')

myRMSE(true_runoff(kk,gi)', ens_AG.meanpriorrunoff(kk,gi)')
myRMSE(true_runoff(kk,gi)', ens_AG.mean_post_runoff(kk,gi)')

%% Plot ensemble of runoff (Ens_LM) with shadederrorbar

% Show median +- IQR
ens_LM.prior_runoff_med = median(ens_LM.prior_runoff_ens,3)';
ens_LM.post_runoff_med = median(ens_LM.post_runoff,3)';
ens_LM.prior_runoff_IQR = iqr(ens_LM.prior_runoff_ens,3)';
ens_LM.post_runoff_IQR = iqr(ens_LM.post_runoff,3)';
ens_LM.prior_runoff_std = std(ens_LM.prior_runoff_ens, [], 3)';
ens_LM.post_runoff_std = std(ens_LM.post_runoff, [], 3)';
% Better yet, show median, 2.5th percentile, and 97.5th percentile
% to reflect UR95 better and avoid making it look like there can be negative values 

ens_LM.prior975 = prctile(ens_LM.prior_runoff_ens, 97.5, 3)';
ens_LM.prior25 = prctile(ens_LM.prior_runoff_ens, 2.5, 3)';
ens_LM.post975 = prctile(ens_LM.post_runoff, 97.5, 3)';
ens_LM.post25 = prctile(ens_LM.post_runoff, 2.5, 3)';

% errbars_prior = [ens_LM.prior975(:,kk), ens_LM.prior25(:,kk)];
% errbars_post = [ens_LM.post975(:,kk), ens_LM.post25(:,kk)];

% errbars_prior = ens_LM.prior_runoff_IQR(:,kk);
% errbars_post = ens_LM.post_runoff_IQR(:,kk);

ens_LM.meanpriorrunoff = mean(ens_LM.prior_runoff_ens,3);
ens_LM.meanpriorrunoff = mean(ens_LM.prior_runoff_ens,3);

figure(1)
subplot(2,2,4)

errbars_prior = 1.96*ens_LM.prior_runoff_std(gi,kk);
errbars_post = 1.96*ens_LM.post_runoff_std(gi,kk);

h1 = shadedErrorBar(tv(gi),ens_LM.meanpriorrunoff(kk,gi),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),ens_LM.mean_post_runoff(kk,gi),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

title(['Runoff at cell ' num2str(kk) ' (Ens LM)'])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
ylim([-6,20])
grid on
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,gi)', ens_LM.meanpriorrunoff(kk,gi)')
myKGE(true_runoff(kk,gi)', ens_LM.mean_post_runoff(kk,gi)')

myNSE(true_runoff(kk,gi)', ens_LM.meanpriorrunoff(kk,gi)')
myNSE(true_runoff(kk,gi)', ens_LM.mean_post_runoff(kk,gi)')

myRMSE(true_runoff(kk,gi)', ens_LM.meanpriorrunoff(kk,gi)')
myRMSE(true_runoff(kk,gi)', ens_LM.mean_post_runoff(kk,gi)')
