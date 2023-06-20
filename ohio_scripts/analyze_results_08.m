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

%% Load outputs from ISR (daily meas)

% load PW13 ISR results
PW13 = load('./ohio_data/results/final/ISR_results_PW13_m240.mat');

% load Y20 ISR results
Y20 = load('./ohio_data/results/final/ISR_results_Y20_m240.mat');

% load ensemble domain ISR results
ENS = load('./ohio_data/results/final/ISR_results_domain_m240.mat');
ENS.mean_post_runoff = mean(ENS.post_runoff,3)';

% need to re-run domain and Y20 no swot with proper gauge errors
% need to fix any NaNs in the Y20 outcomes
% should do runs with much larger L, T

%% Load outputs from ISR (swot meas)

% load PW13 ISR results
PW13 = load('./ohio_data/results/final/ISR_results_PW13_m240_swot_revised.mat');

% load Y20 ISR results
Y20 = load('./ohio_data/results/final/ISR_results_Y20_swot.mat');

% load ensemble domain ISR results
ENS = load('./ohio_data/results/final/ISR_results_domain_m240_swot.mat');
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

tmpa.predQ = state_model_dumb(tmpa_runoff_prior, HH);

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
t = [42, 73, 88, 180];
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

t = [42, 73, 88, 180];
cbnds = [-1,15;-1,6;-1,10;-1,3];

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

% Sort gauges from upstream to downstream (by drainage area)

[gage_area_sorted,ind_inc_DA] = sort(basin.gage_area, 'descend'); % gauges, sorted from ds to us
% Loop over gauges from smallest basin to largest basin
% Populate map with number of competing grid cells

% basin.mask(isnan(basin.mask)) = 0;
% basin.mask = logical(double(basin.mask));
competing_gages = basin.mask;
bmask = basin.mask;
bmask(isnan(bmask))=0;
bmask = logical(bmask);

for mm=1:m
    % for first subbasin,
    current_gage = ind_inc_DA(mm);
    current_mask = flipud(basin.mask_gage(:,:,current_gage));
    ncells_in_subbasin = nansum(current_mask(:));
    
    % get index of current subbasin cells
    
    current_mask_lin = current_mask(:) == 1;
    ind_current_basin = find(current_mask_lin);
    
%     [ii1, ii2] = find(current_mask == 1);
    competing_gages(ind_current_basin) = ncells_in_subbasin;
end

figure
plotraster(basin.lonv, basin.latv, competing_gages, 'Competing grid cells')
% hold on
% plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', 20)
colormap(flipud(parula(50)))

%% Is there a correlation between posterior NSE and #competinggages

for kk=1:n
%     nse_post(kk) = myNSE(nldas_runoff_true(:,kk), PW13.post_runoff_PW13(:,kk));
    nse_post(kk) = myNSE(nldas_runoff_true(:,kk), Y20.post_runoff_Y20(:,kk));
%     nse_post(kk) = myNSE(nldas_runoff_true(:,kk), ENS.mean_post_runoff(:,kk));
    nse_prior(kk) = myNSE(nldas_runoff_true(:,kk), tmpa_runoff_prior(:,kk));
end

nse_post_map = make_map(basin, nse_post);
nse_change_map = make_map(basin, nse_post - nse_prior);

figure
plotraster(basin.lonv, basin.latv, nse_post_map, 'NSE post')
caxis([0,1])

figure
plotraster(basin.lonv, basin.latv, nse_change_map, 'NSE change')

figure
plot(competing_gages(:)/3681, nse_change_map(:), 'k.')
xlabel('# competing')
ylabel('NSE change')

figure
plot(competing_gages(:)/3681, nse_post_map(:), 'k.')
xlabel('# competing')
ylabel('post NSE')

%% Plot map of sub-basins?

% smaller subbasins should have better results
% number the sub-basins, use colorpicker for this

% Might be better done in a mapping program
% Need a raster with different numbers for each subbasin

figure,imagesc(basin.mask_gage(:,:,1))

% trying to figure out how to choose interesting grid cells to plot/explain
% what we see in the results


%% Plot runoff timeseries for selected grid cells

basin.lon = basin.grid_lon(flipud(bmask));
basin.lat = basin.grid_lat(flipud(bmask));

gc = [1410,199, 50, 3000]; % grid cells to plot

% [ii1,ii2] = ind2sub(size(basin.mask), gc(kk));
% 
% % get index in basin_mask of the cell I want to plot
% plotlon = -85.75
% plotlat = 39

gi = gi';
figure
for kk=1:length(gc)
    
    subplot(2,2,kk)
    plot(tv, tmpa_runoff_prior(:, gc(kk)), 'green*', 'linewidth', lw);
    hold on
    plot(tv, PW13.post_runoff_PW13(:, gc(kk)), 'linewidth', lw, 'color', 'blue');
    plot(tv, Y20.post_runoff_Y20(:,gc(kk)), 'linewidth', lw, 'color', 'cyan');
    plot(tv, ENS.mean_post_runoff(:, gc(kk)), 'r--', 'linewidth', lw);
    plot(tv, nldas_runoff_true(:, gc(kk)), 'linewidth', lw, 'color', 'k');
    legend('TMPA (prior)','PW13','Y20','ENS','NLDAS (true)')
    title(['Grid cell' num2str(gc(kk))])
%     xlim([datetime(2009,6,1), datetime(2009,9,1)])
    ylim([-3,30])

    nse_prior = myNSE(nldas_runoff_true(gi, gc(kk)), tmpa_runoff_prior(gi, gc(kk)));
    nse_PW13 = myNSE(nldas_runoff_true(gi, gc(kk)), PW13.post_runoff_PW13(gi, gc(kk)));
    nse_Y20 = myNSE(nldas_runoff_true(gi, gc(kk)), Y20.post_runoff_Y20(gi,gc(kk)));
    nse_ENS = myNSE(nldas_runoff_true(gi, gc(kk)), ENS.mean_post_runoff(gi, gc(kk)));
    text(7, 10, ['Prior NSE: ' num2str(nse_prior)])
    text(7, 15, ['PW13 NSE: ' num2str(nse_PW13)])
    text(7, 20, ['Y20 NSE: ' num2str(nse_Y20)])
    text(7, 25, ['ENS NSE: ' num2str(nse_ENS)])
    set(gca, 'fontsize', fs)
    
end
%%
% Plot a map showing these points on the runoff NSE map
for kk=1:n
    nse_post(kk) = myNSE(nldas_runoff_true(:,kk), Y20.post_runoff_Y20(:,kk));
end
nse_post_map = make_map(basin, nse_post);

figure
plotraster(basin.lonv, basin.latv, nse_post_map, 'nse')
caxis([0,1])
hold on
plot(basin.lon(gc(kk)), basin.lat(gc(kk)), 'r.')



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
% basin_in.mask = flipud(basin_in.mask);

nse_min = 0;
ms = 10;
figure
subplot(2,2,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, tmpa_runoff_prior, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
title('TMPA Prior NSE')
caxis([nse_min,1])

subplot(2,2,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, PW13.post_runoff_PW13, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
title('PW13 Posterior NSE')
caxis([nse_min,1])

subplot(2,2,3)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, Y20.post_runoff_Y20, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
title('Y20 Posterior NSE')
caxis([nse_min,1])

subplot(2,2,4)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin_in, ENS.mean_post_runoff, nldas_runoff_true', gi);
title('ENS Posterior NSE')
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)
caxis([nse_min,1])

%% Reproducing Yang figure 7

% Choose a gage
% Plot its watershed
% Plot its discharge, plus our estimates

g1 = 1;
g2 = 47;
figure
subplot(2,2,1)
plot(tv, true_discharge(:,g1), 'k', 'linewidth', lw)
hold on
plot(tv, PW13.predQ(:,g1), 'blue--', 'linewidth', lw)
plot(tv, Y20.predQ(:,g1), 'red--', 'linewidth', lw)
plot(tv, ENS.predQ(:,g1), 'cyan', 'linewidth', lw)
plot(tv, tmpa.predQ(:,g1), 'green--', 'linewidth', lw)
legend('Truth','PW13','Y20','Ens', 'Prior')
xlabel('Day')
ylabel('Discharge (mm/day)')

subplot(2,2,2)
plotraster(basin.lonv, basin.latv, flipud(basin.mask_gage(:,:,g1)), ['Basin ' num2str(g1)])

subplot(2,2,3)
plot(tv, true_discharge(:,g2), 'k', 'linewidth', lw)
hold on
plot(tv, PW13.predQ(:,g2), 'blue--', 'linewidth', lw)
plot(tv, Y20.predQ(:,g2), 'red--', 'linewidth', lw)
plot(tv, ENS.predQ(:,g2), 'cyan', 'linewidth', lw)
plot(tv, tmpa.predQ(:,g2), 'green--', 'linewidth', lw)
legend('Truth','PW13','Y20','Ens', 'Prior')
legend('Truth')
xlabel('Day')
ylabel('Discharge (mm/day)')

subplot(2,2,4)
plotraster(basin.lonv, basin.latv, flipud(basin.mask_gage(:,:,g2)), ['Basin ' num2str(g2)])

