% Allegheny Y20 vs Ens ISR comparison

%% Load setup

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('/Users/jschap/Documents/MATLAB/')) % optional cbrewer/colorspace (FEX)
addpath(genpath('./src/'))

lw=2;fs=18;

load('./allegheny_data/alleg_nldas_nanfilled.mat')
load('./allegheny_data/alleg_tmpa_nanfilled.mat')
load('./allegheny_data/setup-swot-gage.mat');

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

true_discharge = state_model_dumb(true_runoff', HH);
[nt,m] = size(true_discharge);

%% Load ISR results

% load('./ohio_data/results/final_aug24/m1s2L40T5_prior_ens.mat');
Y20 = load('./allegheny_data/final/Y20_m0s2L40T5_null_daily_add_15.mat');
Y20.postQ = state_model_dumb(Y20.post_runoff_Y20, HH);

% ENS = load('./allegheny_data/ensISR_m0s2L40T5_swot_null_0.mat');
ENS = load('./allegheny_data/final/ens_M3000_m0s2L40T5_daily_add_15.mat');
% ENS.mean_post_runoff = mean(ENS.post_runoff,3)';
ENS.predQ = state_model_dumb(ENS.mean_post_runoff', HH);
% ENS.median_post_runoff = median(ENS.post_runoff,3)'; % median vs. mean is
% similar here
% ENS.predQ = state_model_dumb(ENS.median_post_runoff, HH);

%% Assemble true and prior mean runoff matrices

figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;
tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

totmean = mean(ENS.truth.total_runoff(:));
null_runoff_prior = totmean*ones(size(tmpa_runoff_prior));
yprior = state_model_dumb(null_runoff_prior, HH);
% yprior = state_model_dumb(tmpa_runoff_prior, HH);

%% Figure 1. Runoff snapshot maps.

% Where do neg val occur?
Y20.minval = min(Y20.post_runoff_Y20, [], 2);
ii=find(Y20.minval<0);
sum(ii)
figure,plot(Y20.minval)

sum(Y20.post_runoff_Y20(:)<0) % number of negative values
100*sum(Y20.post_runoff_Y20(:)<0)/(365*3681) % percent neg values
100*sum(Y20.post_runoff_Y20(Y20.post_runoff_Y20<0))/sum(Y20.post_runoff_Y20(Y20.post_runoff_Y20>0)) % percent neg volume

basin.true_runoff = nldas_runoff_true;
basin.prior_runoff = tmpa_runoff_prior;

t = [75,128];
cbnds = [-1,4;-1,10];

% 102, 120, 123, 124 (candidate days)
% for jj=41:50
%     t(2)=ii(jj);
nse = plot_runoff_map_snapshots(t, cbnds, basin, 0, ...
    Y20.post_runoff_Y20, ENS.mean_post_runoff, tmpa, nldas);
% end
    
% adjusting colorbar to highlight negative values
nn = 21; % must be odd
cmap = cbrewer2('seq','YlGnBu',nn); % blues at bottom
colormap([1,0,0;cmap]); % black for negative values

%savefig(gcf, './ohio_data/results/final_aug24/Figures/for_paper/m1s2L40T5_50gage_tmpa_daily_15ro_snaps2.fig')
%saveas(gca, './ohio_data/results/final_aug24/Figures/for_paper/m1s2L40T5_50gage_tmpa_daily_15ro_snaps2.png')

tv(t)

%% Figure 2. Runoff time series at a cell (not using in the current draft)

addpath /Users/jschap/Documents/MATLAB

mean_prior_ens = mean(prior_ensemble,3);
% Find a cell with good improvement in NSE
nse_prior = zeros(n,1);
nse_post =zeros(n,1);
for k=1:n
    nse_post(k) = myNSE(nldas_runoff_true(:,k), ENS.mean_post_runoff(:,k));
    nse_prior(k) = myNSE(nldas_runoff_true(:,k), mean_prior_ens(k,:)');
end
figure,plot(nse_post)

% 1412, 1520 are good ones to show
cell1 = 1412;
cell2 = 1400;

% Use IQR, not stddev

% [mpost, errpost] = shaded_runoff_plot_lognorm(ENS.post_runoff, cell1);
% [mprior, errprior] = shaded_runoff_plot_lognorm(prior_runoff_ens, cell1);
% % showing lognorm dist
% figure
% plot(1:365, nldas_runoff_true(:, cell1), 'linewidth', lw, 'color', 'k');
% hold on
% shadedErrorBar(1:365,mprior,errprior,'lineprops','-green','patchSaturation',0.33)
% shadedErrorBar(1:365,mpost,errpost,'lineprops',{'-ro','MarkerFaceColor','r'})

prior=struct();
prior.ens_cell1 = squeeze(prior_ensemble(cell1,:,:))';
prior.iq_radius_cell1 = iqr(prior.ens_cell1', 2)/2;
prior.median_cell1 = median(prior.ens_cell1,1)';
prior.errbars_cell1 = [prior.iq_radius_cell1, prior.iq_radius_cell1];

ENS.ens_cell1 = squeeze(ENS.post_runoff(cell1,:,:))';
ENS.iq_radius_cell1 = iqr(ENS.ens_cell1', 2)/2;
ENS.median_cell1 = median(ENS.ens_cell1,1)';
ENS.errbars_cell1 = [ENS.iq_radius_cell1, ENS.iq_radius_cell1];

% mean and stddev
figure
s1 = shadedErrorBar(1:365,prior.median_cell1,prior.errbars_cell1,'lineprops','-blue','patchSaturation',0.33);
s1.mainLine.LineWidth = lw;
hold on
s2=shadedErrorBar(1:365,ENS.median_cell1,ENS.errbars_cell1,'lineprops','-red','patchSaturation',0.075);
s2.mainLine.LineWidth = lw;
plot(1:365, nldas_runoff_true(:, cell1), 'linewidth', lw, 'color', 'k');
legend('Prior','Posterior','Truth')
title('Runoff at Cell 100 (-88.4375, 38.0625)')
xlim([1,365])
% ylim([-10,30])
xlabel('DOY')
ylabel('Runoff (mm/day)')
set(gca, 'fontsize', fs)

% subplot(2,1,2)
% shadedErrorBar(1:365,squeeze(prior_ensemble(cell2,:,:))',{@mean,@std},'lineprops','-green','patchSaturation',0.33)
% hold on
% shadedErrorBar(1:365,squeeze(ENS.post_runoff(cell2,:,:))',{@mean,@std},'lineprops','-r','patchSaturation',0.33)
% plot(1:365, nldas_runoff_true(:, cell2), 'linewidth', lw, 'color', 'k');
% ylim([-10,40])
% xlim([1,365])
% xlabel('DOY')
% ylabel('Runoff (mm/day)')
% legend('Prior','Posterior','Truth')
% title('Runoff at Cell 1400 (-85.4375, 37.9375)')
% set(gca, 'fontsize', fs)

%% Figure 3. Discharge estimates (similar to Y20 plot for subbasins)

% Choose a gage
% Plot its watershed
% Plot its discharge, plus our estimates

% Need to find a gauge in a headwater subbasin:

g1 = ENS.cal_ind(1); % 31 is big

ms = 20;

figure
subplot(1,3,[1,2])
plot(tv, true_discharge(:,g1), 'k', 'linewidth', lw)
hold on
plot(tv, Y20.postQ(:,g1), 'red-*', 'linewidth', lw)
plot(tv, ENS.predQ(:,g1), 'cyan', 'linewidth', lw)
plot(tv, yprior(:,g1), 'green-', 'linewidth', lw)
plot(tv, Y20.gage_w_error(:,g1), 'k.', 'markersize', ms)
xlim([tv(60), tv(120)]) % 151 or 243
legend('Truth','Y20','Ens', 'Prior', 'Obs')
title('Additive errors, m0s2L40T5, null prior')
xlabel('Day')
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)

subplot(1,3,3)
plotraster(basin.lonv, basin.latv, flipud(basin.mask_gage(:,:,g1)), ['Basin ' num2str(g1)])
hold on
plot(basin.gage_lon(g1), basin.gage_lat(g1), 'r.', 'markersize', ms)
myKGE(true_discharge(:,g1), yprior(:,g1))
myKGE(true_discharge(:,g1), Y20.postQ(:,g1))
myKGE(true_discharge(:,g1), ENS.predQ(:,g1))

%savefig(gcf, './ohio_data/results/final_aug24/Figures/for_paper/m1s2L40T5_50gage_tmpa_swot_15_QTS_2.fig')
%saveas(gca, './ohio_data/results/final_aug24/Figures/for_paper/m1s2L40T5_50gage_tmpa_swot_15_QTS_2.png')

%% Calculate KGE for calibration/val gauges

kge=zeros(m,3);
for kk=1:m
    kge(kk,1) = myKGE(true_discharge(:,kk), yprior(:,kk));
    kge(kk,2) = myKGE(true_discharge(:,kk), Y20.postQ(:,kk));
    kge(kk,3) = myKGE(true_discharge(:,kk), ENS.predQ(:,kk));
end

median(kge(cal_ind,:),1)
median(kge(val_ind,:),1)

%% Figure 4. Runoff NSE maps

nse_min = -1;
ms = 15;
cal_ind=ENS.cal_ind;

k=size(HH,3)-1;
gi = 1*(k+1):nt-(k+1);

figure(16)
subplot(1,3,1)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, tmpa_runoff_prior, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon(cal_ind), basin.gage_lat(cal_ind), 'r.', 'markersize', ms)
title('TMPA Prior NSE')
caxis([nse_min,1])

figure(16)
subplot(1,3,2)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, Y20.post_runoff_Y20, nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon(cal_ind), basin.gage_lat(cal_ind), 'r.', 'markersize', ms)
title('Y20 Posterior NSE')
caxis([nse_min,1])

figure(16)
subplot(1,3,3)
[nse, kge, rmse, nsemap] = plot_gofmaps(basin, ENS.mean_post_runoff, nldas_runoff_true', gi);
title('ENS Posterior NSE')
hold on 
plot(basin.gage_lon(cal_ind), basin.gage_lat(cal_ind), 'r.', 'markersize', ms)
caxis([nse_min,1])


sum(ENS.mean_post_runoff(:)<0)/numel(ENS.mean_post_runoff(:))
sum(Y20.post_runoff_Y20(:)<0)/numel(Y20.post_runoff_Y20(:))

% nn = 21; % must be odd
% cmap = cbrewer2('seq','YlGnBu',nn); % blues at bottom
% (cmap); % black for negative values

% greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
% redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
% colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
% colormap(flipud(colorMap))

%saveas(gca, './ohio_data/results/final_aug24/Figures/for_paper/nse_maps2.png')
%savefig(gcf, './ohio_data/results/final_aug24/Figures/for_paper/nse_maps2.fig')

%% Figure 5. Runoff NSE change maps

figure, [nse, kge, rmse, nsemap_tmpa] = plot_gofmaps(basin, tmpa_runoff_prior, nldas_runoff_true', gi);
figure, [nse, kge, rmse, nsemap_Y20] = plot_gofmaps(basin, Y20.post_runoff_Y20, nldas_runoff_true', gi);
figure, [nse, kge, rmse, nsemap_ENS] = plot_gofmaps(basin, ENS.mean_post_runoff, nldas_runoff_true', gi);

figure
subplot(2,2,1)
plotraster(basin.lonv, basin.latv, nsemap_Y20-nsemap_tmpa, 'TMPA \DeltaNSE')
hold on 
plot(basin.gage_lon(Y20.cal_ind), basin.gage_lat(Y20.cal_ind), 'r.', 'markersize', ms)
caxis([-5,5])

subplot(2,2,2)
plotraster(basin.lonv, basin.latv, nsemap_ENS-nsemap_tmpa, 'ENS \DeltaNSE')
hold on 
plot(basin.gage_lon(Y20.cal_ind), basin.gage_lat(Y20.cal_ind), 'r.', 'markersize', ms)
caxis([-5,5])

%saveas(gca, './ohio_data/results/final_aug24/Figures/for_paper/deltanse_maps2.png')
%savefig(gcf, './ohio_data/results/final_aug24/Figures/for_paper/deltanse_maps2.fig')
