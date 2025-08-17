% Make figures for the ISR paper (with log EnsISR (not used, though))
%
% 11/10/2023 JRS

clear
cd('/home/jschap/Dropbox/ISR/inverse_streamflow_routing/allegheny_data/results')
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
y20_AG = load('./outlet_only/Cases_1-4/Case1_Y20AG_m0a1L40T5_daily_5_unbiased_M1-500.mat'); % Case 1 AG errors, Y20 prop
y20_LM = load('./outlet_only/Cases_1-4/Case2_Y20LM_m1a1L40T5_daily_5_unbiased_M1-500.mat'); % Case 2 LM errors, Y20 prop
% y20_LM = load('./outlet_only/Cases_1-4/Case3.1_Y20LM_m0a1L40T5_daily_5_unbiased_flat_alpha4_M1-500.mat'); % Case 3.1 LM errors, Y20 flat
% y20_LM = load('./outlet_only/Cases_1-4/Case3.2_Y20LM_m0a1L40T5_daily_5_unbiased_flat_alpha0.25_M1-500.mat'); % Case 3.2 LM errors, Y20 flat
% y20_TMPA = load('./outlet_only/Cases_1-4/Case4_Y20AG_m0a1L40T5_daily_5_tmpa_outlet_only.mat')

% low prior (0.5*nldas) (0.25*tmpa)
% y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_low_prior.mat');

% medium prior (tmpa)
% y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_med_prior.mat');

% high prior (2*nldas) (4*tmpa)
% y20_AG = load('./allegheny_data/results/outlet_only/Y20AG_m0a1L40T5_daily_5_high_prior.mat');

% AG vs. LM, Y20 vs. Ens
% Unbiased prior, 15% error, 10-day meas, 14 gauges (ens of Y20 --> use
% make_figs_ens mfile)
ens_LM = load('./allgage/Ens/ensLM_m1a1L40T5_swot_15.mat');
ens_AG = load('./allgage/Ens/ensAG_m0a1L40T5_swot_15.mat');
% y20_AG = load('./allegheny_data/results/allgage/Y20AG_m0a1L40T5_daily_0_tmpa_flat.mat');
% y20_LM = load('./allegheny_data/results/outlet_only/Y20LM_m1a1L40T5_daily_5_tmpa_flat.mat');

% Effect of proportional runoff variance assumption
% TMPA, daily, 14 gauges
% y20_AG = load('./allegheny_data/results/allgage/Y20AG_m0a1L40T5_daily_0_tmpa_times4_flt.mat');
% ens_AG = load('./allegheny_data/results/allgage/ensAGflat_m0a1L40T5_swot_15.mat');

y20_AG.prior_runoff_ens = permute(y20_AG.prior_runoff_ens, [2,1,3]);
% y20_LM.prior_runoff_ens = permute(y20_LM.prior_runoff_ens, [2,1,3]);

onerealiz = 1;
if onerealiz
% if just using one realization for Y20
    mm=1;
    y20_AG.post_runoff = y20_AG.post_runoff(:,:,mm);
    y20_AG.post_stddev = y20_AG.post_stddev(:,:,mm);
    y20_AG.prior_stddev = y20_AG.prior_stddev(:,:,mm);
    y20_AG.runoff_prior = y20_AG.prior_runoff_ens(:,:,mm);
    y20_LM.post_runoff = y20_LM.post_runoff(:,:,mm);
    y20_LM.post_stddev = y20_LM.post_stddev(:,:,mm);
    y20_LM.prior_stddev = y20_LM.prior_stddev(:,:,mm);
    y20_LM.runoff_prior = y20_LM.prior_runoff_ens(:,:,mm);    
else
    1;
end

% Y20_m0a1L40T5 (former)

% mean prior runoff
ens_AG.meanpriorrunoff = mean(ens_AG.prior_runoff_ens,3)';
ens_AG.mean_post_runoff = mean(ens_AG.post_runoff,3)';
ens_LM.meanpriorrunoff = mean(ens_LM.prior_runoff_ens,3)';
ens_LM.mean_post_runoff = mean(ens_LM.post_runoff,3)';

% renaming some terms
% y20_AG.runoff_prior = y20_AG.tmpa_runoff_prior;
% y20_AG.post_runoff = y20_AG.post_runoff_Y20;
% y20_AG.tmpa_runoff_prior = y20_AG.runoff_prior;

%% Double check error distributions

y20_AG.runoff_errors = y20_AG.prior_runoff_ens - true_runoff';
y20_LM.runoff_errors = y20_LM.prior_runoff_ens./true_runoff';

figure,histogram(y20_AG.runoff_errors(1,1,:)), title('AG')
figure,histogram(y20_LM.runoff_errors(1,1,:)), title('LM')

%% Show where negative values occur (and how cKF can fix them)

tt = 110;
prior_runoff_snapshot = make_map(basin, y20_daily.tmpa_runoff_prior(tt,:));
post_runoff_snapshot = make_map(basin, y20_daily.post_runoff(tt,:));
post_runoff_snapshot = make_map(basin, y20_daily_cKF.post_runoff(tt,:));
true_runoff_snapshot = make_map(basin, y20_daily_cKF.nldas_runoff_true(tt,:));

figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, m1, 'Prior')
caxis([-1,6])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, post_runoff_snapshot, 'Posterior')
caxis([-1,6])
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, true_runoff_snapshot, 'Truth')
caxis([-1,6])
% y20_daily.post_runoff
    
% Show whether the mass balance is preserved

gi=4:120; 
priorQ = state_model_dumb(y20_daily.tmpa_runoff_prior, HH);
postQ =  state_model_dumb(y20_daily_cKF.post_runoff, HH);
trueQ =  state_model_dumb(y20_daily_cKF.nldas_runoff_true, HH);
sum(priorQ(gi,2))
sum(postQ(gi,2))
sum(trueQ(gi,2))

% no need to show this, doesn't make a big difference

%% Figure 1 - runoff NSE maps

[kge, rmse, nse] = plt_ISR_results_overview(basin, y20_AG.runoff_prior', ...
    y20_AG.post_runoff', truth, tv, gi)
title('Y20 A/G')

[kge, rmse, nse] = plt_ISR_results_overview(basin, y20_LM.runoff_prior', ...
    y20_LM.post_runoff', truth, tv, gi)
title('Y20 L/M')

plt_ISR_results_overview(basin, ens_AG.meanpriorrunoff', ...
    ens_AG.mean_post_runoff', truth, tv, gi)
title('Ens A/G')

plt_ISR_results_overview(basin, ens_LM.meanpriorrunoff', ...
    ens_LM.mean_post_runoff', truth, tv, gi)
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

% methodology: (use mean + 1.96sigma instead of current method, for
% comparison with y20_AG results ---> doesn't make a difference
% top5 = ens_AG.meanpriorrunoff + 1.96*ens_AG.prior_runoff_std;
% bottom5 = ens_AG.meanpriorrunoff - 1.96*ens_AG.prior_runoff_std;
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
boxplot([y20_AG.ur95_prior, y20_AG.ur95_post, y20_LM.ur95_prior, y20_LM.ur95_post, ens_AG.ur95_prior, ens_AG.ur95_post, ens_LM.ur95_prior, ens_LM.ur95_post], 'labels', ...
    {'Y20 AG Prior','Y20 AG Post','Y20 LM Prior','Y20 LM Post','Ens AG Prior', 'Ens AG Post', 'Ens LM Prior', 'Ens LM Post'},...
    'notch', 'on')
title('95% uncertainty ratio')

%% Figure 3 - ER95 boxplots

y20_AG.er95 = myER95_gauss(truth.true_runoff, ...
    y20_AG.post_runoff', y20_AG.post_stddev');

y20_AG.prior_er95 = myER95_gauss(truth.true_runoff, ...
    y20_AG.runoff_prior', y20_AG.prior_stddev');

ens_AG.er95 = myER(truth.true_runoff, ens_AG.post_runoff, 95);
ens_LM.er95 = myER(truth.true_runoff, ens_LM.post_runoff, 95);

ens_AG.ens_prior_er95 = myER(truth.true_runoff, ens_AG.prior_runoff_ens, 95);
ens_LM.ens_prior_er95 = myER(truth.true_runoff, ens_LM.prior_runoff_ens, 95);

figure
boxplot([y20_AG.prior_er95, y20_AG.er95, ens_AG.ens_prior_er95, ens_AG.er95, ens_LM.ens_prior_er95, ens_LM.er95], 'labels', ...
    {'Y20 AG Prior','Y20 AG Post','Ens AG Prior', 'Ens AG Post', 'Ens LM Prior', 'Ens LM Post'},...
    'notch', 'on')
title('95% exceedence ratio')

% these don't make sense. why are the prior ER95 values different? 

%% Set up colors for figures

colors = [0,0,255; % dark blue (Ens prior)
           0, 255, 0; % green (Y20 posterior)
           255,0,0;  % red (Ens posterior)
           0,255,255]/255; % cyan (Y20 prior)
       
%% Figure 4 - discharge time series

true_discharge = state_model_dumb(true_runoff',HH);
ens_AG.priorQ = state_model_dumb(ens_AG.meanpriorrunoff, HH);
ens_AG.postQ = state_model_dumb(ens_AG.mean_post_runoff, HH);
ens_LM.priorQ = state_model_dumb(ens_LM.meanpriorrunoff, HH);
ens_LM.postQ = state_model_dumb(ens_LM.mean_post_runoff, HH);

y20_AG.priorQ = state_model_dumb(y20_AG.runoff_prior, HH);
y20_AG.postQ = state_model_dumb(y20_AG.post_runoff, HH);

gagei = 2; % outlet

y20_AG.kge_prior = myKGE(true_discharge(gi,gagei),y20_AG.priorQ(gi,gagei));
ens_AG.kge_prior = myKGE(true_discharge(gi,gagei),ens_AG.priorQ(gi,gagei));
ens_LM.kge_prior = myKGE(true_discharge(gi,gagei),ens_LM.priorQ(gi,gagei));

% shadedErrorBar(1:365, squeeze(prior_runoff_ens(kk1,:,:))', {@mean, @std}, ...
%     'lineprops', '-r');

figure
subplot(3,1,1) % Y20 AG
plot(y20_AG.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
y20_AG.kge_post = plot_discharge_ts(y20_AG.postQ(gi,gagei), true_discharge(gi,gagei));
title('Y20 A/G')
plot(y20_AG.gage_w_error(gi,gagei), 'k.', 'markersize', ms)
legend('Prior','Truth','Posterior','Observations')
set(gca, 'fontsize', fs)

subplot(3,1,2) % ENS AG
plot(ens_AG.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
ens_AG.kge_post = plot_discharge_ts(ens_AG.postQ(gi,gagei), true_discharge(gi,gagei));
title('Ens A/G')
plot(ens_AG.error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', ms)
set(gca, 'fontsize', fs)
legend('off')

subplot(3,1,3) % ENS LM
plot(ens_LM.priorQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
ens_LM.kge_post = plot_discharge_ts(ens_LM.postQ(gi,gagei), true_discharge(gi,gagei));
title('Ens L/M')
plot(ens_LM.error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', ms)
set(gca, 'fontsize', fs)
legend('off')

%% Low vs. high prior analysis

% Are low priors with multiplicative errors problematic?

% gi=35:45;

% calculate NSE for all cells
nse.posterior;

% find cells where TMPA prior is very low/very high
tr = true_runoff(:,gi)';
prior = y20_AG.tmpa_runoff_prior(gi,:);
post = y20_AG.post_runoff(gi,:);

kk=40; % very high prior
kk=101; % very low prior

% focus on timesteps 35:45 becuase there is often a peak in true runoff there

% for kk=1:10
%     figure
%     plot(y20_AG.tmpa_runoff_prior(gi,kk), 'blue', 'linewidth', lw)
%     hold on
%     plot(y20_AG.post_runoff(gi,kk), 'red', 'linewidth', lw)
%     plot(tr(gi,kk), 'k', 'linewidth', lw)
%     legend('Prior','Posterior','Truth')
%     title(['Cell ' num2str(kk)])
% end

% how high or low is the prior compared to the truth?
% use one-sided ER95 to see if truth is contained by the prior
z = 1.645;
figure
plot(prior(:,k) - z*y20_AG.prior_stddev(gi,k))
hold on,plot(tr(:,k),'k')

% how do I know if prior is big? small?
% If prior plus 1.645*sd is smaller than truth, prior is small
% If prior minus 1.645*sd is larger than truth, prior is large
z=1.645;
small_ind = prior + z*y20_AG.prior_stddev(gi,:) < tr;
% large_ind = ~small_ind;
large_ind = prior - z*y20_AG.prior_stddev(gi,:) > tr;
% small_ind = prior < tr;
large_ind = prior > tr;

figure
imagesc(small_ind) % there are a good number of small prior locations
% % but very few places where the prior is large, by this metric
% 
% % how good are the estimates where the prior is small vs. other estimates?
% myNSE(tr(small_ind), y20_AG.runoff_prior(small_ind))
% myNSE(tr(~small_ind), y20_AG.runoff_prior(~small_ind))
% 
% myKGE(tr(small_ind), y20_AG.post_runoff(small_ind))
% myKGE(tr(~small_ind), y20_AG.post_runoff(~small_ind))

%% Revised scatterplots figure

figure
subplot(1,2,1)
plot(tr(~small_ind), prior(~small_ind), '.red')
hold on
plot(tr(small_ind), prior(small_ind), '.blue')
ylabel('Estimate')
title('Prior runoff')
legend('Not small prior','Small prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'black', 'LineStyle', '--');

subplot(1,2,2)
plot(tr(~small_ind), post(~small_ind), '.red')
hold on
plot(tr(small_ind), post(small_ind), '.blue')
ylabel('Estimate')
title('Posterior runoff')
legend('Not small prior','Small prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'black', 'LineStyle', '--');

notsmallnseprior = myNSE(tr(~small_ind), prior(~small_ind))
notsmallnsepost = myNSE(tr(~small_ind), post(~small_ind))
smallnseprior = myNSE(tr(small_ind), prior(small_ind))
smallnsepost = myNSE(tr(small_ind), post(small_ind))


%% Scatterplots figure (draft)

figure

subplot(2,3,1)
plot(tr(small_ind), prior(small_ind), '.')
ylabel('Estimate')
title('Prior runoff (small prior)')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,2)
plot(tr(large_ind), prior(large_ind), '.')
ylabel('Estimate')
title('Prior runoff (large prior)')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,3)
plot(tr(:), prior(:), '.')
ylabel('Estimate')
title('Prior runoff (all priors)')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,4)
plot(tr(small_ind), post(small_ind), '.')
ylabel('Estimate')
title('Posterior runoff, given small prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,5)
plot(tr(large_ind), post(large_ind), '.')
ylabel('Estimate')
title('Posterior runoff, given large prior')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

subplot(2,3,6)
plot(tr(:), post(:), '.')
ylabel('Estimate')
title('Posterior runoff, given all priors')
xlabel('Truth')
axis([0,20,0,20])
line([0, 15], [0, 15], 'Color', 'red', 'LineStyle', '--');

% Color code by size of prior runoff relative to truth
% Can use darker colors for larger difference between prior and truth
% Use same color scheme for prior and posterior scatterplot
% Include NSE values on plots
% Repeat for AGflat method and compare to see effect of proportional errors

% Does this change for the AGflat case? 
% If prior runoff is small, does NSE improve?
% If prior runoff is large, does NSE improve?

% overall
myNSE(tr(:), prior(:))
myNSE(tr(:), post(:))

% not small (similar to overall)
myNSE(tr(~small_ind), prior(~small_ind))
myNSE(tr(~small_ind), post(~small_ind))

myNSE(tr(small_ind), prior(small_ind))
myNSE(tr(small_ind), post(small_ind))

myNSE(tr(large_ind), prior(large_ind))
myNSE(tr(large_ind), post(large_ind))

%% Plot ensemble of runoff (Y20_AG) with shadederrorbar

kk=93; % has negative runoff
kk=44; 
kk=101;
kk=40;

figure
%subplot(1,3,2)

meas_times = find(~isnan(y20_AG.gage_w_error(gi,2))); 

% Show median/mean plus or minus IQR

errbars_prior = 1.96*y20_AG.prior_stddev(gi,kk); % 
errbars_post = 1.96*y20_AG.post_stddev(gi,kk); % 

h1 = shadedErrorBar(tv(gi),y20_AG.runoff_prior(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),y20_AG.post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',y20_AG.runoff_prior(gi,kk))
nse2 = myNSE(true_runoff(kk,gi)',y20_AG.post_runoff(gi,kk))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

% plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Medium prior'])
% title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

% ylim([-20,30])
% ylim([-1,7])
% ylim([-4,11])
grid on
set(gca, 'fontsize', fs)
% shows 95% confidence bounds

%% Plot ensemble of runoff (Y20_AG) with lines

z = 1.96;

figure
h1 = plot(y20_AG.runoff_prior(:,kk), 'b*-', 'linewidth', lw);
hold on
plot(y20_AG.runoff_prior(:,kk) + z*y20_AG.prior_stddev(:,kk), 'b-', 'linewidth', lw)
plot(y20_AG.runoff_prior(:,kk) - z*y20_AG.prior_stddev(:,kk), 'b-', 'linewidth', lw)
h2 = plot(y20_AG.post_runoff(:,kk), 'r*-', 'linewidth', lw);
plot(y20_AG.post_runoff(:,kk) + z*y20_AG.post_stddev(:,kk), 'r-', 'linewidth', lw)
plot(y20_AG.post_runoff(:,kk) - z*y20_AG.post_stddev(:,kk), 'r-', 'linewidth', lw)
h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);
legend([h1,h2,h3], 'Prior','Posterior','Truth')
xlabel('t')
ylabel('runoff (mm/day)')

title('Y20, A/G')
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,:)', y20_AG.runoff_prior(:,kk))
myKGE(true_runoff(kk,:)', y20_AG.post_runoff(:,kk))

myNSE(true_runoff(kk,:)', y20_AG.runoff_prior(:,kk))
myNSE(true_runoff(kk,:)', y20_AG.post_runoff(:,kk))

myRMSE(true_runoff(kk,:)', y20_AG.runoff_prior(:,kk))
myRMSE(true_runoff(kk,:)', y20_AG.post_runoff(:,kk))

%% Plot ensemble of runoff (Y20_LM) with shadederrorbar

kk=40;

figure
%subplot(2,1,1)

meas_times = find(~isnan(y20_LM.gage_w_error(gi,2))); 

% Show median/mean plus or minus IQR

errbars_prior = 1.96*y20_LM.prior_stddev(gi,kk); % 
errbars_post = 1.96*y20_LM.post_stddev(gi,kk); % 

h1 = shadedErrorBar(tv(gi),y20_LM.runoff_prior(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),y20_LM.post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

nse1 = myNSE(true_runoff(kk,gi)',y20_LM.runoff_prior(gi,kk))
nse2 = myNSE(true_runoff(kk,gi)',y20_LM.post_runoff(gi,kk))

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

% plot([meas_times meas_times], [-5,20], 'k--', 'linewidth', 1)

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Runoff (mm/day)')

% ylim([-6,20])
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

figure
subplot(2,1,2)
kk=93;
%kk=44;

errbars_prior = 1.96*ens_AG.prior_runoff_std(gi,kk);
errbars_post = 1.96*ens_AG.post_runoff_std(gi,kk);

h1 = shadedErrorBar(tv(gi),ens_AG.meanpriorrunoff(gi,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv(gi),ens_AG.mean_post_runoff(gi,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3 = plot(tv(gi), true_runoff(kk,gi), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Mean','Posterior Mean',...
    'Truth');

title(['Runoff at cell ' num2str(kk)])
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
ylim([-6,20])
grid on
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,gi)', ens_AG.meanpriorrunoff(gi,kk))
myKGE(true_runoff(kk,gi)', ens_AG.mean_post_runoff(gi,kk))

myNSE(true_runoff(kk,gi)', ens_AG.meanpriorrunoff(gi,kk))
myNSE(true_runoff(kk,gi)', ens_AG.mean_post_runoff(gi,kk))

myRMSE(true_runoff(kk,gi)', ens_AG.meanpriorrunoff(gi,kk))
myRMSE(true_runoff(kk,gi)', ens_AG.mean_post_runoff(gi,kk))

%% GoF statistics

% Change in runoff RMSE for all cells for each method

%% Plot ensemble of runoff (Ens_LM) with shadederrorbar

% Show median +- IQR
ens_LM.prior_runoff_med = median(ens_LM.prior_runoff_ens,3)';
ens_LM.post_runoff_med = median(ens_LM.post_runoff,3)';
ens_LM.prior_runoff_IQR = iqr(ens_LM.prior_runoff_ens,3)';
ens_LM.post_runoff_IQR = iqr(ens_LM.post_runoff,3)';

% Better yet, show median, 2.5th percentile, and 97.5th percentile
% to reflect UR95 better and avoid making it look like there can be negative values 

errbars_prior = ens_LM.prior_runoff_IQR(:,kk);
errbars_post = ens_LM.post_runoff_IQR(:,kk);

figure
h1 = shadedErrorBar(tv,ens_LM.prior_runoff_med(:,kk),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,ens_LM.post_runoff_med(:,kk),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);

legend([h1.mainLine, h2.mainLine,h3], ...
    'Prior Median','Posterior Median',...
    'Truth');

title('Ens runoff estimates at grid cell 93')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)

% Quantify uncertainty shrinkage using IQR -> UR95

%% Plot ensemble of runoff (Ens_LM) with lines

figure
h1 = plot(ens_LM.meanpriorrunoff(:,kk), 'b*-', 'linewidth', lw);
hold on
plot(ens_LM.prior_runoff_2p5(:,kk), 'b-', 'linewidth', lw)
plot(ens_LM.prior_runoff_97p5(:,kk), 'b-', 'linewidth', lw)
h2 = plot(ens_LM.mean_post_runoff(:,kk), 'r*-', 'linewidth', lw);
plot(ens_LM.post_runoff_2p5(:,kk), 'r-', 'linewidth', lw)
plot(ens_LM.post_runoff_97p5(:,kk), 'r-', 'linewidth', lw)
h3 = plot(true_runoff(kk,:), 'k-', 'linewidth', lw);
legend([h1,h2,h3], 'Prior','Posterior','Truth')
xlabel('t')
ylabel('runoff (mm/day)')
title('Ens, L/M')
set(gca, 'fontsize', fs)

myKGE(true_runoff(kk,:)', ens_LM.meanpriorrunoff(:,kk))
myKGE(true_runoff(kk,:)', ens_LM.mean_post_runoff(:,kk))

myNSE(true_runoff(kk,:)', ens_LM.meanpriorrunoff(:,kk))
myNSE(true_runoff(kk,:)', ens_LM.mean_post_runoff(:,kk))

myRMSE(true_runoff(kk,:)', ens_LM.meanpriorrunoff(:,kk))
myRMSE(true_runoff(kk,:)', ens_LM.mean_post_runoff(:,kk))

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

%% Plot discharge with uncertainty (ENS LM)

% It might be worth expressing uncertainty in terms of quantiles instead of
% standard deviation..

% ensAG
[ens_AG.priorQmean, ens_AG.priorQstd, ens_AG.priorQens, ens_AG.priormedianQ,ens_AG.priorq2p5, ens_AG.priorq97p5, ens_AG.priorQiqr] = ...
    calc_ensemble_discharge_moments(ens_AG.prior_runoff_ens, HH);
[ens_AG.postQmean, ens_AG.postQstd, ens_AG.postQens, ens_AG.postmedianQ,ens_AG.postq2p5, ens_AG.postq97p5, ens_AG.postQiqr] = ...
    calc_ensemble_discharge_moments(ens_AG.post_runoff, HH);

% ensLM
[ens_LM.priorQmean, ens_LM.priorQstd, ens_LM.priorQens, ens_LM.priormedianQ,ens_LM.priorq2p5, ens_LM.priorq97p5, ens_LM.priorQiqr] = ...
    calc_ensemble_discharge_moments(ens_LM.prior_runoff_ens, HH);
[ens_LM.postQmean, ens_LM.postQstd, ens_LM.postQens, ens_LM.postmedianQ,ens_LM.postq2p5, ens_LM.postq97p5, ens_LM.postQiqr] = ...
    calc_ensemble_discharge_moments(ens_LM.post_runoff, HH);

%% Comparing distributions of prior and posterior runoff btw Y20 and EnsISR

% Plot PDF for Y20 prior and posterior runoff

runoff_case = 1;
switch runoff_case
    case 1
        tt=110;
        kk=93; % highlighting a time when Y20 posterior runoff<0 and when stddev is reduced significantly
        x = linspace(-50,150);
        ymax = 0.06;
    case 2
        tt=41;
        kk=47; % highlighting a time when Y20 and EnsISR both give similar (accurate) results
        x = linspace(-5,10);
        ymax = 0.6;
end
xmin = min(x);
xmax = max(x);

figure
y = normpdf(x,y20_AG.runoff_prior(tt,kk),y20_AG.prior_stddev(tt,kk));
plot(x,y, 'color', colors(1,:))
hold on
y = normpdf(x,y20_AG.post_runoff(tt,kk),y20_AG.post_stddev(tt,kk));
plot(x,y, 'color', colors(3,:))
plot([true_runoff(kk,tt), true_runoff(kk,tt)], [0, ymax], 'k--', 'linewidth', lw);
title('Y20 runoff estimates')

figure
histogram(ens_AG.prior_runoff_ens(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(1,:))
hold on
plot([true_runoff(kk,tt), true_runoff(kk,tt)], [0, ymax], 'k--', 'linewidth', lw);
histogram(ens_AG.post_runoff(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(3,:))
title('EnsISR (with log transform) runoff estimates')
xlim([xmin,xmax])

% Plot histogram for EnsISR prior and posterior runoff

figure
histogram(ens_LM.prior_runoff_ens(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(1,:))
hold on
plot([true_runoff(kk,tt), true_runoff(kk,tt)], [0, ymax], 'k--', 'linewidth', lw);
histogram(ens_LM.post_runoff(kk,tt,:), 'normalization', 'pdf', 'facecolor', colors(3,:))
title('EnsISR (with log transform) runoff estimates')
xlim([xmin,xmax])

%% Discharge with confidence bounds (Y20)

gg=2;
lw=3;
fs = 22;
z = 1.96;

errbars_prior = zeros(nt,1);
errbars_post = zeros(nt,1);

% errbars_prior = z*y20_AG.prior_stddevQ(:,gg);
% errbars_post = z*y20_AG.post_stddevQ(:,gg);

tvdates = datetime(2009,1,1):datetime(2009,12,31);
tvdates = tvdates(i1:i2);
tv = 1:120;

figure
h1 = shadedErrorBar(tv,y20_AG.priorQ(:,gg),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

h2 = shadedErrorBar(tv,y20_AG.postQ(:,gg),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(ens_AG.priorQmean(:,gg), '-*', 'linewidth', lw, 'color', colors(1,:));
h4=plot(ens_AG.postQmean(:,gg), '-*', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(true_discharge(:,gg), 'k-', 'linewidth', lw);
h6 = plot(ens_AG.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

% retrace the EnsISR posterior median for visibility
plot(tv,y20_AG.postQ(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));

% xlim([15,85])
legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Y20 Prior', 'Y20 Posterior', ...
    'Ens Prior','Ens Posterior',...
    'Truth','Observations');
% this plot shows median +- IQR. For normal distribution, median = mean

title('Discharge estimates at Allegheny Basin outlet')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)


%% Discharge with confidence bounds (EnsISR)

gg=2;
lw=3;
fs = 22;
z = 1.96;

errbars_prior = z*ens_AG.priorQstd(:,gg);
errbars_post = z*ens_AG.postQstd(:,gg);

% errbars_prior = [ens_AG.priorq97p5(:,gg),ens_LM.priorq2p5(:,gg)];
% errbars_post = [ens_LM.postq97p5(:,gg),ens_LM.postq2p5(:,gg)];
% errbars_prior = ens_LM.priorQiqr(:,gg);
% errbars_post = ens_LM.postQiqr(:,gg);

tvdates = datetime(2009,1,1):datetime(2009,12,31);
tvdates = tvdates(i1:i2);
tv = 1:120;

figure
h1 = shadedErrorBar(tv,ens_AG.priorQmean(:,gg),errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

h2 = shadedErrorBar(tv,ens_AG.postQmean(:,gg),errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(y20_AG.priorQ(:,gg), '-*', 'linewidth', lw, 'color', colors(1,:));
h4=plot(y20_AG.postQ(:,gg), '-*', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(true_discharge(:,gg), 'k-', 'linewidth', lw);
h6 = plot(ens_AG.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

% retrace the EnsISR posterior median for visibility
plot(tv,ens_AG.postQmean(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));

% xlim([15,85])
legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth','Observations');
% this plot shows median +- IQR. For normal distribution, median = mean

title('Discharge estimates at Allegheny Basin outlet')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Discharge, with means and standard deviations for EnsISR

figure

h1=plot(ens_LM.priorQmean(:,gg), '-', 'linewidth', lw, 'color', colors(1,:));
hold on
plot(ens_LM.priorQmean(:,gg) + ens_LM.priorQstd(:,gg), '-', 'linewidth', 1, 'color', colors(1,:));
plot(ens_LM.priorQmean(:,gg) - ens_LM.priorQstd(:,gg), '-', 'linewidth', 1, 'color', colors(1,:));

h2=plot(ens_LM.postQmean(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));
plot(ens_LM.postQmean(:,gg) + ens_LM.postQstd(:,gg), '-', 'linewidth', 1, 'color', colors(3,:));
plot(ens_LM.postQmean(:,gg) - ens_LM.postQstd(:,gg), '-', 'linewidth', 1, 'color', colors(3,:));

h5=plot(y20_AG.priorQ(:,gg), '-*', 'linewidth', lw, 'color', colors(1,:));
h6=plot(y20_AG.postQ(:,gg), '-*', 'linewidth', lw, 'color', colors(2,:));

h3=plot(true_discharge(:,gg), 'k-', 'linewidth', lw);
h4 = plot(ens_LM.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

legend([h5,h6,h1,h2,h3,h4],'Y20 Prior','Y20 Posterior','Ens Prior','Ens Posterior','Truth', 'Obs')
xlabel('t')
% xlim([20,100])
ylim([0,1200])
title('Ens LM outlet discharge (15% error)')
ylabel('discharge (mm/day)')
set(gca, 'fontsize', fs)

%% Show uncertainty ratio for discharge

m=14; % number of gauges
for gg=1:m
    
gi = 4:120; % avoids spin up (no estimate)
top5 = ens_AG.postQmean(gi,gg) + z*ens_AG.postQstd(gi,gg);
bottom5 = ens_AG.postQmean(gi,gg) - z*ens_AG.postQstd(gi,gg);
ens_AG.ur95_post(gg) = sum(top5 - bottom5,1)'./sum(true_discharge(gi,gg));

top5 = ens_AG.priorQmean(gi,gg) + z*ens_AG.priorQstd(gi,gg);
bottom5 = ens_AG.priorQmean(gi,gg) - z*ens_AG.priorQstd(gi,gg);
ens_AG.ur95_prior(gg) = sum(top5 - bottom5,1)'./sum(true_discharge(gi,gg));

end

% Change in UR95 shows how much the 95 percent confidence interval shrunk
delta_UR = ens_AG.ur95_prior - ens_AG.ur95_post
% you want the posterior to be smaller than the prior

figure
boxplot([ens_AG.ur95_prior', ens_AG.ur95_post'])

%% Show 95% exceedence ratio for discharge

ens_AG.er95_prior = myER(true_discharge', permute(ens_AG.postQens,[2,1,3]), 95);
ens_AG.er95_post = myER(true_discharge', permute(ens_AG.priorQens,[2,1,3]), 95);

delta_ER = ens_AG.er95_prior - ens_AG.er95_post
% you want the posterior to be smaller than the prior

figure
boxplot([ens_AG.er95_prior, ens_AG.er95_post])
