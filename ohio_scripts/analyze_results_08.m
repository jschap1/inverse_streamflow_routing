% Look at results and generate figures for publication
%
% 6/8/2023 JRS

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

load('./ohio_data/ohio_tmpa_nanfilled.mat')
load('./ohio_data/ohio_nldas_nanfilled.mat')
load('./ohio_data/swot_like_measurements_100m_no_error_revised.mat')

% add distmat
A = load('./ohio_data/setup-1-gage.mat');
basin.distmat = A.basin.distmat;

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

[nt,m] = size(true_discharge);

%% Load outputs from ISR

% load PW13 ISR results
PW13 = load('./ohio_data/ISR_results_PW13_m240.mat');

% load Y20 ISR results
Y20 = load('./ohio_data/ISR_results_Y20_m240.mat');

% load ensemble domain ISR results
% ENS = load('./ohio_data/domain_ISR_results_M5.mat');

%% Calculate discharge at each gage

PW13.predQ = state_model_dumb(PW13.post_runoff_PW13, HH);
Y20.predQ = state_model_dumb(Y20.post_runoff_Y20, HH);
% ENS.predQ = state_model_dumb(ENS.post_runoff, HH);

%% 

%% Compare prior, posterior, truth (basin mean)

mean_prior_runoff = mean(tmpa_runoff_prior,2);
mean_true_runoff = squeeze(nanmean(nanmean(nldas.runoff,1),2));

% mean_post_runoff = mean(post_runoff_PW13,2);
% mean_post_runoff = mean(post_runoff_Y20,2);
% mean_post_runoff = nanmean(real(ENS.mean_posterior_runoff),1);

lw = 2;
fs = 16;

figure
plot(tv, mean_prior_runoff, 'linewidth', lw, 'color', 'green')
hold on
plot(tv, Y20.mean_post_runoff, 'linewidth', lw, 'color', 'blue')
plot(tv, mean(ENS.mean_posterior_runoff,1), 'r--', 'linewidth', lw)
plot(tv, mean_true_runoff, 'linewidth', lw, 'color', 'k')
legend('Prior','Posterior (Y20)','Posterior (Ensemble ISR)', 'True')
set(gca, 'fontsize', fs)

prior_nse = myNSE(mean_true_runoff, mean_prior_runoff);
post_nse = myNSE(mean_true_runoff, mean_post_runoff');

prior_bias = mean(mean_prior_runoff)/mean(mean_true_runoff);
post_bias = mean(mean_post_runoff)/mean(mean_true_runoff);

%% Posterior runoff timeseries figure (PW13/Y20 method)

figure

subplot(2,1,1)
h1 = plot(tv, mean_true_runoff, 'linewidth', lw, 'color', 'black');
hold on
h2 = plot(tv, mean_prior_runoff, 'b', 'linewidth', lw);
h3 = plot(tv, real(mean_post_runoff), 'r', 'linewidth', lw);
plot([tv(k+1), tv(k+1)],[0 10], 'k')
plot([tv(nt - (k+1)), tv(nt-(k+1))],[0 10], 'k')
grid on
% xlim([tv(i1), tv(i2)])
% xlim([datetime(2006,5,13), datetime(2006,5,14)])
% xlim([datetime(2006,5,13), datetime(2006,5,18)])
xlabel('Time')
ylabel('Basin mean runoff (mm/day)')
legend([h1,h2,h3], {'Truth','Prior Mean','Posterior Mean'})
set(gca, 'fontsize', fs)
ylim([-2,14])
title('ISR')

sdprior = mean(ENS.sd_prior_runoff,1);
sdpost = nanmean(ENS.sd_posterior_runoff,1);

subplot(2,1,2)
h1 = plot(tv, mean_true_runoff, 'linewidth', lw, 'color', 'black');
hold on

h2 = plot(tv, mean_prior_runoff', 'b', 'linewidth', lw);
plot(tv, mean_prior_runoff' + sdprior, '--b', 'linewidth', lw)
plot(tv, mean_prior_runoff' - sdprior, '--b', 'linewidth', lw)

h3 = plot(tv, mean_post_runoff, 'r', 'linewidth', lw);
plot(tv, mean_post_runoff + sdpost, 'r--', 'linewidth', lw)
plot(tv, mean_post_runoff - sdpost, 'r--', 'linewidth', lw)

plot([tv(k+1), tv(k+1)],[0 10], 'k')
plot([tv(nt - (k+1)), tv(nt-(k+1))],[0 10], 'k')
grid on

xlabel('Time')
ylabel('Basin mean runoff (mm/day)')
title('Ensemble ISR')
legend([h1,h2,h3], {'Truth','Prior Mean','Posterior Mean'})
set(gca, 'fontsize', fs)
ylim([-2,14])

%% Compare prior, posterior, and truth (maps)

% Something may be wrong with the posterior_runoff_maps...

% TMPA prior, NLDAS truth

% These should be the same

% figure
% subplot(1,2,1)
% plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), 'TMPA 1')
% subplot(1,2,2)
% plotraster(basin.lonv, basin.latv, make_map(basin, tmpa_runoff_prior(t(i), :)), 'TMPA 1')

% post_runoff_w_neg = post_runoff_PW13;
% post_runoff_PW13(post_runoff_PW13<0) = 0;

% post_runoff_Y20_w_neg = post_runoff_Y20;
% post_runoff_Y20(post_runoff_Y20<0) = 0;

figure
% t=30:34;
t=[42, 75, 81, 336];

cmax = 6;

for i=1:4

%     prior_runoff_map = make_map(basin, tmpa_runoff_prior(t(i), :));
%     post_runoff_map = make_map(basin, post_runoff_Y20(t(i),:));
    post_runoff_map_ENS = make_map(basin, ENS.mean_posterior_runoff(:,t(i)));
    
    % Prior
    subplot(4,4,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['Prior (day ' num2str(t(i)) ')'])
%     plotraster(basin.lonv, basin.latv, prior_runoff_map, ['Prior (day ' num2str(t(i)) ')'])
    caxis([0,cmax])

    % Posterior (Y20)
%     subplot(4,4,4+i)
%     plotraster(basin.lonv, basin.latv, post_runoff_map, ['Posterior (day ' num2str(t(i)) ')'])
%     caxis([0,cmax])
    
    % Posterior (Ensemble)
    subplot(4,4,8+i)
    plotraster(basin.lonv, basin.latv, post_runoff_map_ENS, ['Posterior (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
    
    % Truth
    subplot(4,4,12+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['Truth (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
    
end
colormap cool

%% Plot maps, calculate GoF

% gi = 1:365;
gi = (k+2):nt-(k+1);
% gi = 33:40;

figure
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, tmpa_runoff_prior, nldas_runoff_true', gi);
title('Prior NSE')

figure
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, post_runoff_PW13, nldas_runoff_true', gi);
title('Posterior NSE')

%% Checking Kalman update for perfectly correlated case

sigma = 0.83; % for s = 2 case
C = ones(12,12);
H = [0,0,1,0,1,1,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,1,1,0,1;
    0,0,0,0,0,0,0,0,0,0,1,0];
Cxy = sigma^2*C*H';
Cyy = sigma^2*H*C*H';
Cv = 0.15^2*eye(3);
K = Cxy*inv(Cyy+Cv);