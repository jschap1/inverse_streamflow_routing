% Generates the following figures:

% Figure 4.1.1 Explaining Y20 posterior runoff patterns
% Figure 4.1.2 Runoff time series (E1)
% Figure 4.1.3 Absolute error histograms
% Figure 4.1.5 Negative runoff map snapshot
% Figure 4.2.1 Y20 Q vs. Ens for gage case (outlet/upstream gage)
% Figure 4.2.2 Y20 Q vs. Ens for gage case (outlet/upstream gage)
% Figure 4.2.3 Y20 Q vs. Ens for SWOT case (outlet/upstream gage)
% Figure 4.2.4 Y20 Q vs. Ens for SWOT case (outlet/upstream gage)
% Figure 4.2.5 Runoff NSE and dNSE maps
% Figure 4.2.6 Runoff time series (SWOT case, E4)

% What about the following figures?
% Figure 1.1 Example basin setup for ISR
% Figure 3.1 Study area and data --> generate_figure_3.1.m
% Figure 4.1.4 Two-cell negative runoff case --> generate_figure_4.1.4.m
% Figure 4.3.1 NSE change maps, no correlation - localization
% experiments and/or generate_figures_4_3.m
% Figure 4.3.2 NSE change maps, L40T5 - localization
% experiments and/or generate_figures_4_3.m
% Figure 4.3.3 EnISR with localization Q estimates - use this script,
% repurposes Figure 4.2.4 generator

%% Figures for paper (E1-E4)
%
% 12/8/2023 JRS

clear
cd('/home/jschap/Dropbox/ISR/inverse_streamflow_routing/allegheny_data/results')
addpath(genpath('/home/jschap/Dropbox/ISR/inverse_streamflow_routing/src'))
addpath(genpath('/home/jschap/Documents/MATLAB'))
load('../setup/setup-swot-gage.mat');
i1=60; i2=179;
nt=i2-i1+1;
lw=2;
fs=16;
ms=30;

charlbl =  compose("(%s)",('a':'z').'); % labels

truth = struct();
true_runoff = true_runoff(:,i1:i2);
truth.true_runoff = true_runoff;
truth.total_runoff = true_runoff;
basin.true_runoff = true_runoff;

tv = 1:120;
gi = (k+1):nt-(k+1);

% Load experiments
E1 = load('./E1-8/final/E1_Y20AG_m0a1L40T5_daily_5_tmpa_outletonly.mat'); % Y20, gage
E2 = load('./E1-8/final/E2_ensAG_m0a1L40T5_daily_5_tmpa_M500_outletonly.mat'); % Ens, gage
% E2 = load('./E1-8/final/E4_ensAG_m0a1L40T5_swot_15_tmpa_M500_local.mat'); % E7
E3 = load('./E1-8/final/E3_Y20AG_m0a1L40T5_swot_15_tmpa.mat'); % Y20, SWOT
E4 = load('./E1-8/final/E4_ensAG_m0a1L40T5_swot_15_tmpa_M500.mat'); % Ens, SWOT

% mean prior runoff (Ens)
E2.meanpriorrunoff = mean(E2.prior_runoff_ens,3)';
E2.mean_post_runoff = mean(E2.post_runoff,3)';
E4.meanpriorrunoff = mean(E4.prior_runoff_ens,3)';
E4.mean_post_runoff = mean(E4.post_runoff,3)';

colors = [0,0,255; % dark blue (Ens prior)
           0, 255, 0; % green (Y20 posterior)
           255,0,0;  % red (Ens posterior)
           0,255,255]/255; % cyan (Y20 prior)

%% Figure 4.1.5

% Show where negative values occur (and how cKF can fix them)

tt = 110;
prior_runoff_snapshot = make_map(basin, E1.runoff_prior(tt,:)');
post_runoff_snapshot = make_map(basin, E1.post_runoff(tt,:)');
true_runoff_snapshot = make_map(basin, true_runoff(:,tt));

clim = [-2,4];
figure

subplot(1,3,1)
plotraster(basin.lonv, basin.latv, prior_runoff_snapshot, 'Prior runoff (mm)')
caxis(clim)
a = gca; % or whatever you use to access your axis handle
a.CLim = max(abs(a.CLim)) * [-1 1];
text(0.025,0.05,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left

subplot(1,3,2)
plotraster(basin.lonv, basin.latv, post_runoff_snapshot, 'Posterior runoff (mm)')
caxis(clim)
a = gca; % or whatever you use to access your axis handle
a.CLim = max(abs(a.CLim)) * [-1 1];
text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left

subplot(1,3,3)
plotraster(basin.lonv, basin.latv, true_runoff_snapshot, 'True runoff (mm)')
caxis(clim)
a = gca; % or whatever you use to access your axis handle
a.CLim = max(abs(a.CLim)) * [-1 1];
text(0.025,0.05,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left

greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
colormap(bluewhitered)
% colormap(flipud(cool));

% y20_daily.post_runoff
    
% see if there is a threshold that we can use to predict negative runoff

% for tt=1:120
%     rat(tt) = prctile(E1.runoff_prior(tt,:),95)/prctile(E1.runoff_prior(tt,:),50);
%     minv(tt) = min(E1.post_runoff(tt,:));
% end
% 
% figure,plot(minv, rat, 'o') % the relationship should be more straightforward in the no-corr case...

%% Runoff NSE maps

basin.true_runoff = true_runoff;

[kge, rmse, nse] = plt_ISR_results_overview(basin, E1.runoff_prior', ...
    E1.post_runoff', truth, tv, gi)
title('E1 (gage, Y20)')

[kge, rmse, nse] = plt_ISR_results_overview(basin, E2.meanpriorrunoff', ...
    E2.mean_post_runoff', truth, tv, gi)
title('E2 (gage, Ens)')

plt_ISR_results_overview(basin, E3.runoff_prior', ...
    E3.post_runoff', truth, tv, gi)
title('E3 (SWOT, Y20)')

plt_ISR_results_overview(basin, E4.meanpriorrunoff', ...
    E4.mean_post_runoff', truth, tv, gi)
title('E4 (SWOT, Ens)')
       
%% could the Ens cases have higher spatial correlation than the Y20 cases?

E1.runoff_errors = E1.runoff_prior - true_runoff';
E2.runoff_errors = permute((E2.prior_runoff_ens - true_runoff), [2,1,3]);

mean(E1.runoff_errors(:))
mean(E2.runoff_errors(:))

std(E1.runoff_errors(:))
std(E2.runoff_errors(:)) % E2 runoff errors are twice as large as E1 runoff errors
% try regenerating the ensemble of runoff errors for the Ens cases
% is the difference that Ens errors have std proportional to true runoff
% vs. runoff errors? This could play a role. 

mean(true_runoff(:))^2
mean(E1.tmpa_runoff_prior(:))^2


%% Check that NSE is being calculated correctly

kk=40;
figure
plot(E3.runoff_prior(:,kk))
hold on
plot(E4.meanpriorrunoff(:,kk))
plot(true_runoff(kk,:), 'k')
legend('Y20','Ens', 'true')
title('prior')

figure
plot(E3.post_runoff(:,kk))
hold on
plot(E4.mean_post_runoff(:,kk))
plot(true_runoff(kk,:), 'k')
legend('Y20','Ens', 'true')
title('posterior')

for kk=1:233
    E1.nse(kk) = myNSE(true_runoff(kk,gi)', E1.post_runoff(gi,kk));
    E2.nse(kk) = myNSE(true_runoff(kk,gi)', E2.post_runoff(gi,kk));
    E3.nse(kk) = myNSE(true_runoff(kk,gi)', E3.post_runoff(gi,kk));
    E4.nse(kk) = myNSE(true_runoff(kk,gi)', E4.mean_post_runoff(gi,kk));
end

%% Calculate discharge

B = load('../../allegheny_data/setup/setup-233-gage.mat');
basin = B.basin;
gage = B.gage;
A = B.HH;
clearvars B

% truth
true_discharge = state_model_dumb(true_runoff',A);

% y20 cases
E1.priorQ = state_model_dumb(E1.runoff_prior, A);
E3.priorQ = state_model_dumb(E3.runoff_prior, A);
E1.postQ = state_model_dumb(E1.post_runoff, A);
E3.postQ = state_model_dumb(E3.post_runoff, A);

% ens cases
[E2.meanQprior, E2.stdQprior] = calc_ensemble_discharge_moments(E2.prior_runoff_ens, A);
[E2.meanQpost, E2.stdQpost] = calc_ensemble_discharge_moments(E2.post_runoff, A);
[E4.meanQprior, E4.stdQprior] = calc_ensemble_discharge_moments(E4.prior_runoff_ens, A);
[E4.meanQpost, E4.stdQpost] = calc_ensemble_discharge_moments(E4.post_runoff, A);

% Calculate KGE at river gauges

% Get grid cells on the river channel above a threshold flowacc
figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up
basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;
facc_linear = reshape(flipud(basin.flow_accum), length(basin_mask_linear), 1);
facc = facc_linear(logical(basin_mask_linear),:)';
facc(isnan(facc)) = 0;
thres = 5000e6;
rivind = facc(:)>thres; % rivind == 29, 146 for outlet and upstream gage (ii=4,26)
gageinds = find(rivind);

basin_area = max(basin.gage_area)/1000^2; % basin area = 33751 km2
f = (basin_area/n)*1000/86400; % cubic meters per second conversion factor

for ii=1:length(gageinds)

    gagei = gageinds(ii); % outlet (cell 29), upstream (cell 146)

    % note: NSE is not unit-dependent...
    E1.NSE_prior(ii) = myNSE(f*true_discharge(gi,gagei),f*E1.priorQ(gi,gagei));
    E3.NSE_prior(ii) = myNSE(f*true_discharge(gi,gagei),f*E3.priorQ(gi,gagei));
    E2.NSE_prior(ii) = myNSE(f*true_discharge(gi,gagei),f*E2.meanQprior(gi,gagei));
    E4.NSE_prior(ii) = myNSE(f*true_discharge(gi,gagei),f*E4.meanQprior(gi,gagei));
    
    E1.NSE_post(ii) = myNSE(f*true_discharge(gi,gagei),f*E1.postQ(gi,gagei));
    E3.NSE_post(ii) = myNSE(f*true_discharge(gi,gagei),f*E3.postQ(gi,gagei));
    E2.NSE_post(ii) = myNSE(f*true_discharge(gi,gagei),f*E2.meanQpost(gi,gagei));
    E4.NSE_post(ii) = myNSE(f*true_discharge(gi,gagei),f*E4.meanQpost(gi,gagei));

    % Calculate discharge UR95 for prior vs. posterior to quantify
    % uncertainty reduction
   
end

% Might be handy to have a table of prior/posterior NSE at gages
% from each scenario - let's do it

mean(E1.NSE_prior)
mean(E2.NSE_prior)
mean(E3.NSE_prior)
mean(E4.NSE_prior)
mean(E1.NSE_post)
mean(E2.NSE_post)
mean(E3.NSE_post)
mean(E4.NSE_post)
% we should also calculate this for the localized case (we will case it E5)

%% Plot true and predicted discharge at each gauge

% 2x1 figures showing discharge at the each gauge. One figure for gage
% case (E1,2), one figure for SWOT case (E3,4)
% OR
% A single 2x2 figure, with labeled panels

lw=3;
fs = 20;
z = 1.96;
ms = 25;
tvdates = datetime(2009,1,1):datetime(2009,12,31);
tvdates = tvdates(i1:i2);
tv = 1:120;

% gg=29; % outlet
% gg=146; % upstream gauge

%% Figure 4.2.1/4.2.2 - Y20 vs. Ens for gage case (outlet/upstream gage)

gg=146; % outlet gage or upstream gage (29, 146)

% Convert mm/day to cms for figures
basin_area = max(basin.gage_area)/1000^2; % basin area = 33751 km2
f = (basin_area/n)*1000/86400; % cubic meters per second conversion factor

errbars_prior = z*E2.stdQprior(:,gg);
errbars_post = z*E2.stdQpost(:,gg);

figure

subplot(2,2,1) % Y20, full
plot(f*E1.priorQ(:,gg), '-', 'linewidth', lw, 'color', colors(1,:));
hold on
plot(f*E1.postQ(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));
plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
line([20,20],[-100,3000], 'color', 'k', 'linewidth', 1)
line([60,60],[-100,3000], 'color', 'k', 'linewidth', 1)
legend('Prior','Posterior','Truth');
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('E1 (Y20) upstream gage')
% title('E1 (Y20) outlet gage')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)

ylim([-100,700])
% ylim([0,3000])

subplot(2,2,3) % Ens, full
h1 = shadedErrorBar(tv,f*E2.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E2.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
% legend([h1.mainLine, h2.mainLine,h5], ...
%     'Ens Prior', 'Ens Posterior', ...
%     'Truth');
line([20,20],[-100,3000], 'color', 'k', 'linewidth', 1)
line([60,60],[-100,3000], 'color', 'k', 'linewidth', 1)
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('E2 (Ens) upstream gage')
% title('E2 (Ens) outlet gage')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
ylim([-100,700])
% ylim([0,3000])

subplot(2,2,2) % Y20, zoomed
plot(f*E1.priorQ(:,gg), '-', 'linewidth', lw, 'color', colors(1,:));
hold on
plot(f*E1.postQ(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));
plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
% legend('Y20 Prior','Y20 Posterior','Truth');
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('E1 (Y20, zoomed in)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
xlim([20,60])
% ylim([0,1050])
ylim([-30,300])

subplot(2,2,4) % Ens, zoomed
h1 = shadedErrorBar(tv,f*E2.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E2.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
% legend([h1.mainLine, h2.mainLine,h5], ...
%     'Ens Prior', 'Ens Posterior', ...
%     'Truth');
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('E2 (Ens, zoomed in)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
xlim([20,60])
% ylim([0,1050])
ylim([-30,300])

%% Figure 4.2.3/4.2.4 - Y20 vs. Ens for SWOT case (outlet/upstream gage)

gg=29; % outlet gage or upstream gage (29, 146)

% Convert mm/day to cms for figures
basin_area = max(basin.gage_area)/1000^2; % basin area = 33751 km2
f = (basin_area/n)*1000/86400; % cubic meters per second conversion factor

errbars_prior = z*E4.stdQprior(:,gg);
errbars_post = z*E4.stdQpost(:,gg);

meas_times = find(~isnan(E3.gage_w_error(:, 2)));

figure

ylims = [0,3000];
% ylims = [-100,700];

subplot(2,2,1) % Y20, full
h1=plot(f*E3.priorQ(:,gg), '-', 'linewidth', lw, 'color', colors(1,:));
hold on
h2=plot(f*E3.postQ(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));
h3=plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
line([20,20],[-100,3000], 'color', 'k', 'linewidth', 1)
line([60,60],[-100,3000], 'color', 'k', 'linewidth', 1)
line([meas_times, meas_times], [ylims(1),ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
legend([h1,h2,h3],'Prior','Posterior','Truth');
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
% title('E3 (Y20) upstream gage')
title('E3 (Y20) outlet gage')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
ylim(ylims)

subplot(2,2,3) % Ens, full
h1 = shadedErrorBar(tv,f*E4.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E4.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
% legend([h1.mainLine, h2.mainLine,h5], ...
%     'Ens Prior', 'Ens Posterior', ...
%     'Truth');
% plot(f*E4.error_corrupted_discharge_meas(:,2), 'k.','markersize', ms)
line([20,20],[-100,3000], 'color', 'k', 'linewidth', 1)
line([60,60],[-100,3000], 'color', 'k', 'linewidth', 1)
line([meas_times, meas_times], [ylims(1),ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
% title('E4 (Ens) upstream gage')
title('E4 (Ens) outlet gage')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
% ylim([-100,700])
ylim(ylims)

subplot(2,2,2) % Y20, zoomed
plot(f*E3.priorQ(:,gg), '-', 'linewidth', lw, 'color', colors(1,:));
hold on
plot(f*E3.postQ(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));
plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
% legend('Y20 Prior','Y20 Posterior','Truth');
line([meas_times, meas_times], [ylims(1),ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('E3 (Y20, zoomed in)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
xlim([20,60])
ylim([0,1050])
% ylim([-30,300])

subplot(2,2,4) % Ens, zoomed
h1 = shadedErrorBar(tv,f*E4.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E4.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
% legend([h1.mainLine, h2.mainLine,h5], ...
%     'Ens Prior', 'Ens Posterior', ...
%     'Truth');
% plot(f*E4.error_corrupted_discharge_meas(:,2), 'k.','markersize', ms)
line([meas_times, meas_times], [ylims(1),ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('E4 (Ens, zoomed in)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
xlim([20,60])
ylim([0,1050])
% ylim([-30,300])

%% Figure 4.3.2 (revised) - effect of localization on discharge estimates

gg=29; % outlet gage or upstream gage (29, 146)

% Convert mm/day to cms for figures
basin_area = max(basin.gage_area)/1000^2; % basin area = 33751 km2
f = (basin_area/n)*1000/86400; % cubic meters per second conversion factor

meas_times = find(~isnan(E3.gage_w_error(:, 2)));

figure

gg=29;
errbars_prior = z*E2.stdQprior(:,gg);
errbars_post = z*E2.stdQpost(:,gg);
subplot(2,2,1) % Ens, zoomed
h1 = shadedErrorBar(tv,f*E2.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E2.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
h3=plot(f*E2.error_corrupted_discharge_meas(:,2), 'k.','markersize', ms);
line([meas_times, meas_times], [-100,ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('Outlet gage EnISR-local (E7)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
% xlim([20,60])
ylim([0,3200])

gg=146;
errbars_prior = z*E2.stdQprior(:,gg);
errbars_post = z*E2.stdQpost(:,gg);
subplot(2,2,2) % Ens, zoomed
h1 = shadedErrorBar(tv,f*E2.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E2.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
line([meas_times, meas_times], [-100,ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('Upstream gage EnISR-local (E7)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
% xlim([20,60])
ylim([-100,700])

gg=29;
errbars_prior = z*E4.stdQprior(:,gg);
errbars_post = z*E4.stdQpost(:,gg);
subplot(2,2,3) % Ens, zoomed
h1 = shadedErrorBar(tv,f*E4.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E4.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
h3=plot(f*E4.error_corrupted_discharge_meas(:,2), 'k.','markersize', ms);
h4 = line([meas_times, meas_times], [-100,ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('Outlet gage EnISR (E4)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
legend([h1.mainLine, h2.mainLine,h5, h3, h4(1)], ...
    'Prior', 'Posterior', ...
    'Truth', 'Measurements', 'Meas. times','fontsize', 14);
% xlim([20,60])
ylim([0,3200])

gg=146;
errbars_prior = z*E4.stdQprior(:,gg);
errbars_post = z*E4.stdQpost(:,gg);
subplot(2,2,4) % Ens, zoomed
h1 = shadedErrorBar(tv,f*E4.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});
hold on
h2 = shadedErrorBar(tv,f*E4.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);
h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
line([meas_times, meas_times], [-100,ylims(2)], 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
text(0.025,0.95,charlbl{4},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
title('Upstream gage EnISR (E4)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)
% xlim([20,60])
ylim([-100,700])



%% Panel 1b)

gg=29; % outlet

subplot(1,2,2)

errbars_prior = z*E2.stdQprior(:,gg);
errbars_post = z*E2.stdQpost(:,gg);

h1 = shadedErrorBar(tv,f*E2.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

h2 = shadedErrorBar(tv,f*E2.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(f*E1.priorQ(:,gg), '--', 'linewidth', lw, 'color', colors(1,:));
h4=plot(f*E1.postQ(:,gg), '--', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);
h6 = plot(f*E2.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

% retrace the EnsISR posterior median for visibility
% plot(tv,E3.meanQpost(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));

% xlim([15,85])
legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth','Observations');
% this plot shows median +- IQR. For normal distribution, median = mean

text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs)
title('Discharge at outlet (E1, E2)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)

%% Panel 2a)

gg=146; % outlet

figure
subplot(1,2,1)

errbars_prior = z*E4.stdQprior(:,gg);
errbars_post = z*E4.stdQpost(:,gg);

h1 = shadedErrorBar(tv,f*E4.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

myNSE(f*true_discharge(:,gg), f*E4.meanQprior(:,gg))
myNSE(f*true_discharge(:,gg), f*E4.meanQpost(:,gg))
nanmean(f*E4.stdQprior(:,gg))
nanmean(f*E4.stdQpost(:,gg))

h2 = shadedErrorBar(tv,f*E4.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(f*E3.priorQ(:,gg), '--', 'linewidth', lw, 'color', colors(1,:));
h4=plot(f*E3.postQ(:,gg), '--', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);

if gg==29 % we only have measurements at the outlet, not at the upstream gage (146)
    h6 = plot(f*E4.error_corrupted_discharge_meas(:,2), 'k.', 'markersize', ms);
    legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth','Observations');
else
    legend([h1.mainLine, h2.mainLine,h3,h4,h5], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth');
end

text(0.025,0.05,charlbl{1},'Units','normalized','FontSize',fs)
title('Discharge at upstream gauge (E3, E4)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)

%% Panel 2b)

gg=29; % outlet

subplot(1,2,2)

errbars_prior = z*E4.stdQprior(:,gg);
errbars_post = z*E4.stdQpost(:,gg);

h1 = shadedErrorBar(tv,f*E4.meanQprior(:,gg),f*errbars_prior, 'lineprops', ...
    {'color', colors(1,:),'linewidth', lw});

hold on

h2 = shadedErrorBar(tv,f*E4.meanQpost(:,gg),f*errbars_post, 'lineprops', ...
    {'color', colors(3,:),'linewidth', lw}, ...
    'transparent', true);

h3=plot(f*E3.priorQ(:,gg), '--', 'linewidth', lw, 'color', colors(1,:));
h4=plot(f*E3.postQ(:,gg), '--', 'linewidth', lw, 'color', colors(2,:));

h5 = plot(f*true_discharge(:,gg), 'k-', 'linewidth', lw);

if gg==29 % we only have measurements at the outlet, not at the upstream gage (146)
    h6 = plot(f*E4.error_corrupted_discharge_meas(:,2), 'k.', 'markersize', ms);
    legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth','Observations');
else
    legend([h1.mainLine, h2.mainLine,h3,h4,h5], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth');
end

% h6 = plot(f*E4.error_corrupted_discharge_meas(:,gg), 'k.', 'markersize', ms);

% retrace the EnsISR posterior median for visibility
% plot(tv,E3.meanQpost(:,gg), '-', 'linewidth', lw, 'color', colors(3,:));

% xlim([15,85])
legend([h1.mainLine, h2.mainLine,h3,h4,h5,h6], ...
    'Ens Prior', 'Ens Posterior', ...
    'Y20 Prior','Y20 Posterior',...
    'Truth','Observations');
% this plot shows median +- IQR. For normal distribution, median = mean

text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs)
title('Discharge at outlet (E3, E4)')
xlabel('Day (3/2009 - 6/2009)')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', fs)

%% If necessary

% E2.post_runoff = permute(E2.post_runoff, [2,1,3]);
% E4.post_runoff = permute(E4.post_runoff, [2,1,3]);
% 
% % Get indices of grid cells of interest
% % basin.gage_lat
% % lons = [-78.5, -79.25, -79.5];
% % lats = [41.88, 41, 40.38];
% % cell_ind = zeros(3,1);
% % for ii=1:3
% %     cell_ind(ii) = find(abs(basin.gage_lon - lons(ii))<res/2 & abs(basin.gage_lat - lats(ii))<res/2);
% % end
% cell_ind = [208, 100, 59]; % found: cells 1, 2, and 3 of interest
% 
% % choose cells such that prior is bigger, smaller, and about equal to NLDAS
% 
% [v1, c1] = max(mean(E1.tmpa_runoff_prior,1)' - mean(true_runoff,2)) % prior is bigger
% [v2, c2] = min(mean(E1.tmpa_runoff_prior,1)' - mean(true_runoff,2)) % prior is smaller
% 
% [v3, c3] = abs(mean(E1.tmpa_runoff_prior,1)' - mean(true_runoff,2)<0.1) % prior is about the same
% cell_ind = [c1,c2,71];
% 
% % for which cells is there a lot of bias, or not much bias?
% bias = E1.tmpa_runoff_prior - true_runoff';
% up95 = E1.tmpa_runoff_prior + 1.96*E1.prior_stddev;
% lo95 = E1.tmpa_runoff_prior - 1.96*E1.prior_stddev;
% er95 = sum()/n
% 
% figure
% plot(up95(:,1))
% hold on
% plot(lo95(:,1))
% plot(true_runoff(1,:))
% legend('up95','lo95','truth')
% 
% er = myER95_gauss(true_runoff, E1.tmpa_runoff_prior', E1.prior_stddev');
% figure,
% plot(er, E2.nse_post - E2.nse_prior, '.')
% 
% cell_ind = [40, 211, 100]; % final
% cell_ind = [211, 100]; % final final

%% Figure 4.2.5 Runoff maps

% Result:
% median NSE: 0.057192
% median NSE: 0.33517
% median NSE: 0.33715
% median NSE: 0.057192
% median NSE: 0.1196
% median NSE: 0.10801
% median NSE: 0.067606 % E4 local

cell_ind = [23, 155, 211];

res=1/8;

% add labels for cells 1-3
B = load('../setup/setup-swot-gage.mat');
fs2 = 14;

figure
% Prior, first row
subplot(2,3,1)

[nse.prior, kge.prior, rmse.prior, E1.nsemap_prior] = plot_gofmaps(basin, E1.runoff_prior, truth.total_runoff, gi);
hold on
% h1 = plot(basin.gage_lon(gageinds)-res/2, basin.gage_lat(gageinds)-res/2, '.g', 'markersize', ms);
% h1 = plot(basin.gage_lon(gageinds) - res/2, basin.gage_lat(gageinds) - res/2, 'r.', 'markersize', ms);
h2 = plot(basin.gage_lon(29) - res/2, basin.gage_lat(29) - res/2, 'r.', 'markersize', ms); % outlet gage
h3 = plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms); % upstream gage
h4 = plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms); % cells 1-3
text(0.025,0.05,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
legend([h2,h3,h4], 'Measurement locations', 'Upstream gage', 'Cells of interest')
xlabel('Lon')
ylabel('Lat')
title('Prior NSE')

% % Y20, daily measurements (E1)
% subplot(2,5,2)
% [nse.prior, kge.prior, rmse.prior, E1.nsemap_post] = plot_gofmaps(basin, E1.post_runoff, truth.total_runoff, gi);
% hold on
% plot(basin.gage_lon(gageinds)-res/2, basin.gage_lat(gageinds)-res/2, '.g', 'markersize', ms)
% plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.r', 'markersize', ms)
% text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'r', 'fontsize', fs2)
% text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'r', 'fontsize', fs2)
% xlabel('Lon')
% ylabel('Lat')
% title('Posterior, Y20, Gauge')

subplot(2,3,2) % E2
[nse.prior, kge.prior, rmse.prior, E2.nsemap_post] = plot_gofmaps(basin, E2.mean_post_runoff, truth.total_runoff, gi);
hold on
plot(basin.gage_lon(29) - res/2, basin.gage_lat(29) - res/2, 'r.', 'markersize', ms); % outlet gage
plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms); % upstream gage
plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms) % cells 1-3
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
xlabel('Lon')
ylabel('Lat')
title('E2 NSE (EnISR, Gage)')

% Prior, second row 
% subplot(2,5,6)
% [nse.prior, kge.prior, rmse.prior, nsemap] = plot_gofmaps(basin, E3.runoff_prior, truth.total_runoff, gi);
% hold on
% plot(B.basin.gage_lon-res/2, B.basin.gage_lat-res/2, '.g', 'markersize', ms)
% plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.r', 'markersize', ms)
% text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'r', 'fontsize', fs2)
% text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'r', 'fontsize', fs2)
% % text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'r', 'fontsize', fs2)
% title('Prior NSE')
% xlabel('Lon')
% ylabel('Lat')

% subplot(2,3,3) % E3
% [nse.prior, kge.prior, rmse.prior, nsemap] = plot_gofmaps(basin, E3.post_runoff, truth.total_runoff, gi);
% hold on
% plot(B.basin.gage_lon-res/2, B.basin.gage_lat-res/2, '.g', 'markersize', ms)
% plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.r', 'markersize', ms)
% text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'r', 'fontsize', fs2)
% text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'r', 'fontsize', fs2)
% xlabel('Lon')
% ylabel('Lat')
% title('Posterior, Y20, SWOT')

subplot(2,3,3) % E4
[nse.prior, kge.prior, rmse.prior, nsemap] = plot_gofmaps(basin, E4.mean_post_runoff, truth.total_runoff, gi);
hold on
plot(B.basin.gage_lon-res/2, B.basin.gage_lat-res/2, '.r', 'markersize', ms) % swot gages
plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms); % upstream gage
plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms) % cells
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
text(0.025,0.05,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
xlabel('Lon')
ylabel('Lat')
title('E4 NSE (EnISR, SWOT)')

% plot(basin.gage_lon(29) - res/2, basin.gage_lat(29) - res/2, 'r.', 'markersize', ms); % outlet gage
% plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms); % upstream gage
% plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms) % cells 1-3
% text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
% text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
% text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)

% Runoff NSE change map (E2, gage)
subplot(2,3,5) % E2
[nse.prior, kge.prior, rmse.prior, E2.nsemap_prior] = plot_gofmaps(basin, E2.meanpriorrunoff, truth.total_runoff, gi);
[nse.prior, kge.prior, rmse.prior, E2.nsemap_post] = plot_gofmaps(basin, E2.mean_post_runoff, truth.total_runoff, gi);
plotraster(basin.lonv, basin.latv, E2.nsemap_post - E2.nsemap_prior, '\Delta NSE');
caxis([-3,3])
colormap(gca,bluewhitered)
hold on
plot(basin.gage_lon(29) - res/2, basin.gage_lat(29) - res/2, 'r.', 'markersize', ms); % outlet gage
plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms); % upstream gage
plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms) % cells 1-3
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
text(0.025,0.05,charlbl{4},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
xlabel('Lon')
ylabel('Lat')
title('E2 \DeltaNSE, (EnISR, Gage)')

% Runoff NSE change map (E4, SWOT)
subplot(2,3,6) % E4
[nse.prior, kge.prior, rmse.prior, E4.nsemap_prior] = plot_gofmaps(basin, E4.meanpriorrunoff, truth.total_runoff, gi);
[nse.prior, kge.prior, rmse.prior, E4.nsemap_post] = plot_gofmaps(basin, E4.mean_post_runoff, truth.total_runoff, gi);
plotraster(basin.lonv, basin.latv, E4.nsemap_post - E4.nsemap_prior, '\Delta NSE');
caxis([-3,3])
colormap(gca,bluewhitered)
hold on
plot(B.basin.gage_lon-res/2, B.basin.gage_lat-res/2, '.r', 'markersize', ms) % swot gages
plot(basin.gage_lon(146) - res/2, basin.gage_lat(146) - res/2, 'green.', 'markersize', ms); % upstream gage
plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms) % cells 1-3
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
text(0.025,0.05,charlbl{5},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
xlabel('Lon')
ylabel('Lat')
title('E4 \DeltaNSE, (EnISR, SWOT)')

%% Effect of proportional runoff assumption

[nse.prior, kge.prior, rmse.prior, E1.nsemap_post] = plot_gofmaps(basin, E1.post_runoff, truth.total_runoff, gi);
E1.true_runoff = true_runoff';

figure
plotraster(basin.lonv, basin.latv, E1.nsemap_post - E1.nsemap_prior, 'Post - Prior NSE (E1)')
colormap(bluewhitered)

E1.priorbias = E1.runoff_prior(:) - E1.true_runoff(:);
E1.postbias = E1.post_runoff(:) - E1.true_runoff(:);

E1.bias = E1.tmpa_runoff_prior - true_runoff';
figure
histogram(E1.bias)

% find where the true runoff is 2sigma greater than the prior runoff
smallind = E1.true_runoff > E1.runoff_prior + 2*1*E1.runoff_prior; 
largeind = E1.true_runoff < E1.runoff_prior; 

%lowp = prctile(E1.bias(:),5);
%highp = prctile(E1.bias(:),95);

% smallind = E1.bias(:)<lowp;
% largeind = E1.bias(:)>highp;

E1.true_runoff = true_runoff';

figure
plot(E1.true_runoff(:), E1.post_runoff(:), '.')
hold on
h1=plot(E1.true_runoff(smallind), E1.post_runoff(smallind), 'r.');
h2=plot(E1.true_runoff(largeind), E1.post_runoff(largeind), 'g.');
line([0, 15], [0, 15], 'Color', 'k', 'LineStyle', '--');
legend([h2,h1], 'High prior','Low prior')
ylabel('Posterior runoff (mm/day)')
title('E1')
xlabel('Truth (mm/day)')
axis([0,20,0,20])
set(gca, 'fontsize', fs)

other_ind = ~(largeind+smallind);
figure
h2=plot(E1.priorbias(largeind), E1.postbias(largeind), 'cyan.', 'markersize', 10);
hold on
plot(E1.priorbias(other_ind), E1.postbias(other_ind), 'k.', 'markersize', 10);
h1=plot(E1.priorbias(smallind), E1.postbias(smallind), 'r.', 'markersize', 10);
line([-15, 20], [-15, 20], 'Color', 'k', 'LineStyle', '-');
legend([h2,h1], 'High prior','Low prior')
ylabel('Posterior error (mm/day)')
text(15,-3,'Improved estimate', 'fontsize', fs)
text(-10,10,'Degraded estimate', 'fontsize', fs)
title('Posterior vs. prior runoff error')
grid on
xlabel('Prior error (mm/day)')
axis([-15,41,-15,19])
set(gca, 'fontsize', fs)

%% Not used
% 
% % a better question might be, did our NSE get better or worse?
% % given that a cell had too high or too low a prior
% 
% load('../../allegheny_data/setup/alleg_nldas_nanfilled.mat')
% load('../../allegheny_data/setup/alleg_tmpa_nanfilled.mat')
% nldas.runoff  = nldas.runoff(:, :, i1:i2);
% tmpa.runoff  = tmpa.runoff(:, :, i1:i2);
% 
% biasmap = mean(tmpa.runoff - nldas.runoff,3);
% figure
% plotraster(basin.lonv, basin.latv, biasmap, 'TMPA - NLDAS')
% lowp = prctile(biasmap(:), 10);
% highp = prctile(biasmap(:), 90);
% smallind = biasmap<lowp;
% largeind = biasmap>highp;
% 
% figure
% plot(biasmap(:), E1.nsemap_post(:) - E1.nsemap_prior(:), '.')
% hold on
% plot(biasmap(smallind), E1.nsemap_post(smallind) - E1.nsemap_prior(smallind), '.r')
% plot(biasmap(largeind), E1.nsemap_post(largeind) - E1.nsemap_prior(largeind), '.g')
% line([-2, 2], [0, 0], 'Color', 'k', 'LineStyle', '-');
% xlabel('Bias of prior (mm/day)')
% ylabel('Change in NSE')
% 
% biasmap_post = make_map(basin, mean(E1.post_runoff - E1.true_runoff,1)');
% figure
% plotraster(basin.lonv, basin.latv, biasmap_post, 'Posterior bias')
% % Posterior bias looks just like prior bias...
% 
% figure
% plot(biasmap(:), biasmap_post(:), '.')
% hold on
% plot(biasmap(smallind), biasmap_post(smallind), '.r')
% plot(biasmap(largeind), biasmap_post(largeind), '.g')
% line([-2, 2], [-2, 2], 'Color', 'k', 'LineStyle', '-');
% xlabel('Bias of prior (mm/day)')
% ylabel('Bias of posterior (mm/day)')

%% Figure 4.1.1 Explaining Y20 posterior runoff patterns

cell_ind = [23, 155, 211];

figure
subplot(1,2,1)
load('../../allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('../../allegheny_data/setup/alleg_tmpa_nanfilled.mat')
nldas.runoff  = nldas.runoff(:, :, i1:i2);
tmpa.runoff  = tmpa.runoff(:, :, i1:i2);
plotraster(basin.lonv, basin.latv, mean(tmpa.runoff - nldas.runoff,3), 'Prior - true runoff (mm/day)')
% plotraster(basin.lonv, basin.latv, mean(nldas.runoff - tmpa.runoff,3), 'Truth - Prior (mm/day)')
hold on
% plot(basin.gage_lon(gageinds)-res/2, basin.gage_lat(gageinds)-res/2, '.g', 'markersize', ms)
plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms)
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
text(0.025,0.05,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left

E1.rmse_prior = zeros(n,1);
E1.rmse_post = zeros(n,1);
for i=1:n
    E1.rmse_prior(i) = myRMSE(E1.true_runoff(gi,i), E1.runoff_prior(gi,i));
    E1.rmse_post(i) = myRMSE(E1.true_runoff(gi,i), E1.post_runoff(gi,i));
end

E1.rmsemap_prior = make_map(basin, E1.rmse_prior);
E1.rmsemap_post = make_map(basin, E1.rmse_post);

subplot(1,2,2)
plotraster(basin.lonv, basin.latv, E1.nsemap_post - E1.nsemap_prior, '\DeltaNSE')
caxis([-3,3])
hold on
% plot(basin.gage_lon(gageinds)-res/2, basin.gage_lat(gageinds)-res/2, '.g', 'markersize', ms)
plot(basin.gage_lon(cell_ind)-res/2, basin.gage_lat(cell_ind)-res/2, '.k', 'markersize', ms)
text(basin.gage_lon(cell_ind(1)), basin.gage_lat(cell_ind(1)), 'Cell 1', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(2)), basin.gage_lat(cell_ind(2)), 'Cell 2', 'color', 'k', 'fontsize', fs2)
text(basin.gage_lon(cell_ind(3)), basin.gage_lat(cell_ind(3)), 'Cell 3', 'color', 'k', 'fontsize', fs2)
text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
colormap(bluewhitered)

E1.nse_prior = zeros(n,1);
E1.nse_post = zeros(n,1);
for i=1:n
    E1.nse_prior(i) = myNSE(E1.true_runoff(gi,i), E1.runoff_prior(gi,i));
    E1.nse_post(i) = myNSE(E1.true_runoff(gi,i), E1.post_runoff(gi,i));
end

%subplot(2,3,3) % maybe replace this with subplots of runoff at individual cells

% if we plot RMSE maps, then RMSE can be thought similarly to absolute
% error for a single location and time
% might be better for the figure

% truth-prior vs. improvement in NSE

% in general, dark red becomes lighter, as runoff that is low bias either
% stays the same if too low, or gets corrected (show both)
% blue stays blue, as when prior is too large, posterior overcorrects (show
% with inidividal runoff figure)
% scatterplot shows that small gets no change and large can get better or
% worse, show overall stats that posterior runoff error is lower than prior
% runoff error (compute NSE overall)
% very dark red becomes white, since there is no update

myNSE(E1.true_runoff(:), E1.runoff_prior(:))
myNSE(E1.true_runoff(:), E1.post_runoff(:))

% Calculate tabulated values of NSE change for each case in Table 2 (paper)
E1.dnse = E1.nse_post - E1.nse_prior;

ind_larger = mean(E1.true_runoff,1) < mean(E1.runoff_prior,1);
diffvals_tmp = mean(E1.true_runoff,1) - mean(E1.runoff_prior,1);
diffvals_tmp(ind_larger)
med50 = median(diffvals_tmp(ind_larger));
ind_much_larger = diffvals_tmp<med50;
ind_somewhat_larger = ind_larger - ind_much_larger;
ind_much_smaller = mean(E1.true_runoff,1) > 2*mean(E1.runoff_prior,1);
ind_others = ones(1,n) - ind_much_larger - ind_much_smaller - ind_somewhat_larger;
sum(ind_others)+sum(ind_much_smaller) + sum(ind_somewhat_larger) + sum(ind_much_larger) % sums to 233

mean(E1.dnse(logical(ind_somewhat_larger)))
mean(E1.dnse(logical(ind_much_larger)))

ind_others = ones(1,n) - ind_larger - ind_much_smaller;

mean(E1.dnse(logical(ind_much_smaller)))
mean(E1.dnse(logical(ind_larger)))
mean(E1.dnse(logical(ind_others)))

mean_true_runoff = mean(E1.true_runoff,1);
mean_prior_runoff = mean(E1.runoff_prior,1);

mean_prior_runoff(logical(ind_somewhat_larger)) - mean_true_runoff(logical(ind_somewhat_larger))
mean_prior_runoff(logical(ind_much_larger)) - mean_true_runoff(logical(ind_much_larger))

%% Figure 4.1.3 Absolute error histograms

E1.abspriorbias = abs(E1.runoff_prior(:) - E1.true_runoff(:));
E1.abspostbias = abs(E1.post_runoff(:) - E1.true_runoff(:));

figure

subplot(3,3,[1 8])
other_ind = ~(largeind+smallind);
h1=plot(E1.abspriorbias(largeind), E1.abspostbias(largeind), 'r.', 'markersize', 10);
hold on
h2=plot(E1.abspriorbias(other_ind), E1.abspostbias(other_ind), 'cyan.', 'markersize', 10);
h3=plot(E1.abspriorbias(smallind), E1.abspostbias(smallind), 'b.', 'markersize', 10);
line([-15, 20], [-15, 20], 'Color', 'k', 'LineStyle', '-');
legend([h3,h2,h1], 'Prior << truth','Prior < truth','Prior > truth', 'location', 'best')
ylabel('|Posterior - Truth| (mm/day)')
text(7,3,'Improved estimate', 'fontsize', 14)
text(2,8,'Degraded estimate', 'fontsize', 14)
title('Absolute error of prior vs. posterior runoff')
grid on
xlabel('|Prior - Truth| (mm/day)')
axis([0,12,0,12])
set(gca, 'fontsize', fs)

E1.abserr_prior = abs(E1.runoff_prior(:) - E1.true_runoff(:));
E1.abserr_post = abs(E1.post_runoff(:) - E1.true_runoff(:));
E1.dabserr = E1.abserr_post - E1.abserr_prior;

xlims = [-2,2];
% figure
subplot(3,3,3)
histogram(E1.dabserr(smallind), 'facecolor', 'b')
hold on
text(-1.7,250,'Decrease in error')
text(0.3,250,'Increase in error')
ylabel('Frequency')
title('Prior << Truth')
xlabel('\Deltaabsolute error (mm/day)')
xlim(xlims)
set(gca, 'fontsize', fs)

subplot(3,3,6)
histogram(E1.dabserr(other_ind))
ylabel('Frequency')
title('Prior < Truth')
xlabel('\Deltaabsolute error (mm/day)')
xlim(xlims)
set(gca, 'fontsize', fs)

subplot(3,3,9)
histogram(E1.dabserr(largeind), 'facecolor', 'r')
ylabel('Frequency')
title('Prior > Truth')
xlabel('\Deltaabsolute error (mm/day)')
xlim(xlims)
set(gca, 'fontsize', fs)

mean(E1.dabserr(largeind))
mean(E1.dabserr(smallind))
mean(E1.dabserr(other_ind))

%% Runoff time series

% Make this a 3x2 figure showing runoff ts at 3 grid cells
% Use something more subtle to show measurement times for 10-day cases

% Top row is Y20 
% Bottom row is Ens

% E2.post_runoff_std = E2.post_runoff_std';

%% Figure 4.1.2 Runoff time series (E1)

% Define the values of kk
nspec = length(cell_ind); % number of special cells to look at
kk_values = cell_ind(1:nspec);

% Create a figure
figure

% Loop over each kk value

for i = 1:length(cell_ind)
    kk = kk_values(i);
    
%     % Determine subplot position
    switch i
        case 1
            subplot(1,3,1)
        case 2
            subplot(1,3,2)
        case 3
            subplot(1,3,3)
        case 4
            subplot(2,3,6)
    end
    
    priorm = mean(E1.runoff_prior(gi,kk)); % this is the mean prior runoff
    postm = mean(E1.post_runoff(gi,kk)); % this is the mean prior runoff
    truem = mean(true_runoff(kk, gi)); % this is the mean true runoff
    line([0,120], [priorm,priorm], 'color','blue', 'linewidth', 2)
    hold on
    line([0,120], [postm,postm], 'color', 'r', 'linewidth', 2)
    line([0,120], [truem,truem], 'color', 'k', 'linewidth', 2)
        
%     % Figure 2a) and 2c)
    meas_times = find(~isnan(E1.gage_w_error(gi, 2)));
    errbars_prior = 1.96 * E1.prior_stddev(gi, kk);
    errbars_post = 1.96 * E1.post_stddev(gi, kk);

    h1 = shadedErrorBar(tv(gi), E1.runoff_prior(gi, kk), errbars_prior, 'lineprops', ...
        {'color', colors(1, :), 'linewidth', lw});
    hold on
    h2 = shadedErrorBar(tv(gi), E1.post_runoff(gi, kk), errbars_post, 'lineprops', ...
        {'color', colors(3, :), 'linewidth', lw}, ...
        'transparent', true);

    nse1 = myNSE(true_runoff(kk, gi)', E1.runoff_prior(gi, kk));
    nse2 = myNSE(true_runoff(kk, gi)', E1.post_runoff(gi, kk));

    h3 = plot(tv(gi), true_runoff(kk, gi), 'k-', 'linewidth', lw);
    
%     % Add a label
    switch i
        case 1
            title('Cell 1, Prior > Truth')
            text(0.025,0.05,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
        case 2
            title('Cell 2, Prior < Truth')
            text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
        case 3
            title('Cell 3, Prior << Truth')
            text(0.025,0.05,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
        case 4
            text(0.025,0.05,charlbl{6},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
    end
    
    if i==3
        legend([h1.mainLine, h2.mainLine, h3], ...
            'Prior Mean', 'Posterior Mean', ...
            'Truth');
    end
    %title(['Cell ' num2str(i)])
    xlabel('Day (3/2009 - 6/2009)')
    ylabel('Runoff (mm/day)')
    ylim([-5,14])
    xlim([0,120])
    set(gca, 'fontsize', fs)

    % Figure 2b) and 2d)
%     subplot(1, nspec, i);
% 
%     if i==2
%         legend([h1.mainLine, h2.mainLine, h3], ...
%             'Prior Mean', 'Posterior Mean', ...
%             'Truth');
%     end
% 
%     E2.prior_runoff_std = std(E2.prior_runoff_ens, [], 3)';
%     E2.post_runoff_std = std(E2.post_runoff, [], 3)';
% 
%     errbars_prior = 1.96 * E2.prior_runoff_std(gi, kk);
%     errbars_post = 1.96 * E2.post_runoff_std(gi, kk);
% 
%     h1 = shadedErrorBar(tv(gi), E2.meanpriorrunoff(gi, kk), errbars_prior, 'lineprops', ...
%         {'color', colors(1, :), 'linewidth', lw});
%     hold on
%     h2 = shadedErrorBar(tv(gi), E2.mean_post_runoff(gi, kk), errbars_post, 'lineprops', ...
%         {'color', colors(3, :), 'linewidth', lw}, ...
%         'transparent', true);
% 
%     h3 = plot(tv(gi), true_runoff(kk, gi), 'k-', 'linewidth', lw);
% 
%     myNSE(true_runoff(kk, gi)', E2.meanpriorrunoff(gi, kk))
%     myNSE(true_runoff(kk, gi)', E2.mean_post_runoff(gi, kk))
%     
%     title(['Cell ' num2str(i) ' runoff (E1)' ])
%     
%     xlabel('Day (3/2009 - 6/2009)')
%     ylabel('Runoff (mm/day)')
%     set(gca, 'fontsize', fs)
%     ylim([-10,20])
end

%% Figure 4.2.6 Runoff time series (SWOT case, E4)

% Define the values of kk
kk_values = cell_ind(1:nspec);

% Create a figure
figure;

% Loop over each kk value
for i = 1:nspec
    kk = kk_values(i);
    
%     Determine subplot position based on row and column
%     subplot(2, 3, i);
%     
%     Figure 2a) and 2c)
    meas_times = find(~isnan(E3.gage_w_error(gi, 2)));
%     errbars_prior = 1.96 * E3.prior_stddev(gi, kk);
%     errbars_post = 1.96 * E3.post_stddev(gi, kk);
% 
%     h1 = shadedErrorBar(tv(gi), E3.runoff_prior(gi, kk), errbars_prior, 'lineprops', ...
%         {'color', colors(1, :), 'linewidth', lw});
%     hold on
%     
%     Plot meas times
%     ylim = get(gca, 'YLim');
%     line([meas_times, meas_times], ylim, 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
%     
%     h2 = shadedErrorBar(tv(gi), E3.post_runoff(gi, kk), errbars_post, 'lineprops', ...
%         {'color', colors(3, :), 'linewidth', lw}, ...
%         'transparent', true);
% 
%     nsE3 = myNSE(true_runoff(kk, gi)', E3.runoff_prior(gi, kk));
%     nsE4 = myNSE(true_runoff(kk, gi)', E3.post_runoff(gi, kk));
% 
%     h3 = plot(tv(gi), true_runoff(kk, gi), 'k-', 'linewidth', lw);
% 
%     if i==2
%         legend([h1.mainLine, h2.mainLine, h3], ...
%             'Prior Mean', 'Posterior Mean', ...
%             'Truth');
%     end
%         
%     title(['Y20 runoff (E3, Cell ' num2str(i) ')'])
%     xlabel('Day (3/2009 - 6/2009)')
%     ylabel('Runoff (mm/day)')
% 
%     set(gca, 'fontsize', fs)

    % Figure 2b) and 2d)
    subplot(1, nspec, i);

    E4.prior_runoff_std = std(E4.prior_runoff_ens, [], 3)';
    E4.post_runoff_std = std(E4.post_runoff, [], 3)';

    errbars_prior = 1.96 * E4.prior_runoff_std(gi, kk);
    errbars_post = 1.96 * E4.post_runoff_std(gi, kk);

    h1 = shadedErrorBar(tv(gi), E4.meanpriorrunoff(gi, kk), errbars_prior, 'lineprops', ...
        {'color', colors(1, :), 'linewidth', lw});
    hold on
        
    h2 = shadedErrorBar(tv(gi), E4.mean_post_runoff(gi, kk), errbars_post, 'lineprops', ...
        {'color', colors(3, :), 'linewidth', lw}, ...
        'transparent', true);

    h3 = plot(tv(gi), true_runoff(kk, gi), 'k-', 'linewidth', lw);

    % Plot meas times
    %ylim = get(gca, 'YLim');
    line([meas_times, meas_times], ylim, 'Color', [0.9,0.9,0.9], 'LineStyle', '-', 'LineWidth', 1);
    
    if i==3
        legend([h1.mainLine, h2.mainLine, h3], ...
            'Prior Mean', 'Posterior Mean', ...
            'Truth');
    end

        switch i
        case 1
            title('Cell 1, Prior > Truth')
            text(0.025,0.05,charlbl{1},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
        case 2
            title('Cell 2, Prior < Truth')
            text(0.025,0.05,charlbl{2},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
        case 3
            title('Cell 3, Prior << Truth')
            text(0.025,0.05,charlbl{3},'Units','normalized','FontSize',fs) % 0.05 for bottom-left
        end
    
%     title(['Ens runoff (E4, Cell ' num2str(i) ')'])
    xlabel('Day (3/2009 - 6/2009)')
    ylabel('Runoff (mm/day)')
    ylim([-5,14])
    xlim([0,120])
    set(gca, 'fontsize', fs)
end

%% Plot maps of prior and true runoff

trm = make_map(basin, mean(true_runoff,2));
prm = make_map(basin, mean(E1.tmpa_runoff_prior,1)');

figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, prm, 'Prior runoff')
caxis([0.6,2.3])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, trm, 'True runoff')
caxis([0.6,2.3])
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, prm - trm, 'prior - true')

% Is  there a relationship between the posterior NSE and the size of the
% prior?

for kk=1:233
    E2.nse_prior(kk) = myNSE(true_runoff(kk, gi)', E2.meanpriorrunoff(gi, kk));
    E2.nse_post(kk) = myNSE(true_runoff(kk, gi)', E2.mean_post_runoff(gi, kk));
    E2.rmse_prior(kk) = myRMSE(true_runoff(kk, gi)', E2.meanpriorrunoff(gi, kk));
    E2.rmse_post(kk) = myRMSE(true_runoff(kk, gi)', E2.mean_post_runoff(gi, kk));
end

%%

figure
plot(mean(E1.tmpa_runoff_prior,1)' - mean(true_runoff,2), E2.nse_post - E2.nse_prior, '.' ,'markersize', 12)
hold on
plot([0 0], ylim, 'k')
plot(xlim, [0 0], 'k')
title('\Delta NSE vs. (prior - truth)')
xlabel('Prior - Truth (mm/day)')
ylabel('Posterior NSE')

%%

figure
subplot(1,2,1)
plot(mean(E1.tmpa_runoff_prior,1)' - mean(true_runoff,2), E2.nse_prior, '.' ,'markersize', 5)
hold on
axis([-1,1,-5,2])
plot([0 0], ylim, 'k')
plot(xlim, [0 0], 'k')
title('Prior NSE vs. (prior - truth)')
xlabel('Prior - Truth (mm/day)')
ylabel('Posterior NSE')

subplot(1,2,2)
plot(mean(E1.tmpa_runoff_prior,1)' - mean(true_runoff,2), E2.nse_post, '.' ,'markersize', 5)
hold on
axis([-1,1,-5,2])
plot([0 0], ylim, 'k')
plot(xlim, [0 0], 'k')
title('Post NSE vs. (prior - truth)')
xlabel('Prior - Truth (mm/day)')
ylabel('Posterior NSE')

%%
figure
plot(E2.nse_prior, E2.nse_post, '.' ,'markersize', 5)
hold on
plot([0 0], ylim, 'k')
plot(xlim, [0 0], 'k')
xlabel('Prior')
ylabel('Posterior')

%%
figure
plot(mean(true_runoff,2), E2.nse_post, '.' ,'markersize', 5)
hold on
plot([0 0], ylim, 'k')
plot(xlim, [0 0], 'k')
title('NSE vs. prior')
xlabel('Prior (mm/day)')
ylabel('Posterior NSE')

figure
plot(mean(E1.tmpa_runoff_prior,1)', E2.nse_post, '.' ,'markersize', 5)
hold on
plot([0 0], ylim, 'k')
plot(xlim, [0 0], 'k')
title('NSE vs. truth')
xlabel('Truth (mm/day)')
ylabel('Posterior NSE')
