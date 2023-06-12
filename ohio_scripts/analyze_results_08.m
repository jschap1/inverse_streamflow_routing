% Look at results and generate figures for publication
%
% 6/8/2023 JRS

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('/Users/jschap/Documents/MATLAB/')) % optional cbrewer/colorspace (FEX)
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
ENS = load('./ohio_data/ISR_results_domain_m240.mat');
ENS.mean_post_runoff = mean(ENS.post_runoff,3)';

%% Calculate discharge at each gage

PW13.predQ = state_model_dumb(PW13.post_runoff_PW13, HH);
Y20.predQ = state_model_dumb(Y20.post_runoff_Y20, HH);
ENS.predQ = state_model_dumb(ENS.mean_post_runoff, HH);

%% Assemble true and prior mean runoff matrices

basin.mask = flipud(basin.mask);
figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;
tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

%% Compare prior, posterior, truth (basin mean)

mean_prior_runoff = mean(tmpa_runoff_prior,2);
mean_true_runoff = squeeze(nanmean(nanmean(nldas.runoff,1),2));

% mean_post_runoff = mean(post_runoff_PW13,2);
% mean_post_runoff = mean(post_runoff_Y20,2);
% mean_post_runoff = nanmean(real(ENS.mean_posterior_runoff),1);

lw = 2;
fs = 16;

figure
h1 = plot(tv, mean_prior_runoff, 'linewidth', lw, 'color', 'green');
hold on
h2 = plot(tv, mean(PW13.post_runoff_PW13,2), 'linewidth', lw, 'color', 'blue');
h3 = plot(tv, mean(Y20.post_runoff_Y20,2), 'linewidth', lw, 'color', 'cyan');
h4 = plot(tv, mean(ENS.mean_post_runoff,2), 'r--', 'linewidth', lw);
h5 = plot(tv, mean_true_runoff, 'linewidth', lw, 'color', 'k');
set(gca, 'fontsize', fs)

% add vertical lines at key times
t = [42, 73, 81, 336];
ymax = 14;
plot([tv(t(1)) tv(t(1))],[0 ymax], 'k-')
plot([tv(t(2)) tv(t(2))],[0 ymax], 'k-')
plot([tv(t(3)) tv(t(3))],[0 ymax], 'k-')
plot([tv(t(4)) tv(t(4))],[0 ymax], 'k-')

legend([h1, h2, h3, h4, h5], 'Prior','Posterior (PW13)', 'Posterior (Y20)','Posterior (Ensemble ISR)', 'True')

prior_nse = myNSE(mean_true_runoff, mean_prior_runoff);
PW13.post_nse = myNSE(mean_true_runoff, mean(PW13.post_runoff_PW13,2));
Y20.post_nse = myNSE(mean_true_runoff, mean(Y20.post_runoff_Y20,2));
% ENS.post_nse = myNSE(mean_true_runoff, mean(ENS.post_runoff,2));

prior_bias = mean(mean_prior_runoff)/mean(mean_true_runoff);
PW13.post_bias = mean(mean(PW13.post_runoff_PW13,2))/mean(mean_true_runoff);
Y20.post_bias = mean(mean(Y20.post_runoff_Y20,2))/mean(mean_true_runoff);

%% Compare prior, posterior, and true runoff (maps)

% Where do neg val occur?
PW13.minval = min(PW13.post_runoff_PW13, [], 2);
Y20.minval = min(Y20.post_runoff_Y20, [], 2);

t=[42, 73, 88, 336];
cbnds = [-1,20;-1,6;-1,6;-1,10];

plot_runoff_map_snapshots(t, cbnds, basin, PW13.post_runoff_PW13, ...
    Y20.post_runoff_Y20, ENS.mean_post_runoff, tmpa, nldas);

% adjusting colorbar to highlight negative values
nn = 21; % must be odd
cmap = cbrewer2('seq','YlGnBu',nn); % blues at bottom
colormap([1,0,0;cmap]); % black for negative values

%% Plot map of how many downstream gauges each cell has?

basin.ds_gages = flipud(sum(basin.mask_gage,3));

figure
plotraster(basin.lonv, basin.latv, basin.ds_gages, 'ds gages')

%% Plot map of how many other cells this cell has to share a gage with?

%% Plot map of sub-basins?

% smaller subbasins should have better results
% number the sub-basins, use colorpicker for this

% Might be better done in a mapping program
% Need a raster with different numbers for each subbasin

figure,imagesc(basin.mask_gage(:,:,1))

% trying to figure out how to choose interesting grid cells to plot/explain
% what we see in the results


%% Plot runoff timeseries for selected grid cells

gc = []; % grid cells to plot

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

%% Plot maps, calculate GoF

% gi = 1:365;
gi = (k+2):nt-(k+1);
% gi = 33:40;

basin_in = basin;
basin_in.mask = flipud(basin_in.mask);

nse_min = 0;
ms = 10;
figure
subplot(1,4,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, tmpa_runoff_prior, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
title('TMPA Prior NSE')
caxis([nse_min,1])

subplot(1,4,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, PW13.post_runoff_PW13, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
title('PW13 Posterior NSE')
caxis([nse_min,1])

subplot(1,4,3)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, Y20.post_runoff_Y20, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
title('Y20 Posterior NSE')
caxis([nse_min,1])

subplot(1,4,4)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, ENS.mean_post_runoff, nldas_runoff_true', gi);
title('ENS Posterior NSE')
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
caxis([nse_min,1])


