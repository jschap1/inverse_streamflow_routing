% Generate Figures 4.3 for the ISR paper
%
% 12/8/2023 JRS
% Generates Figures 4.3.1/2/3
%
% Material comes from allegheny_scripts/localization_experiments.m

clear, clc, close all
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
cd /home/jschap/Dropbox/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

fs=18;lw=2;
load('./allegheny_data/setup/setup-swot-gage.mat');
aa=load('./allegheny_data/setup/setup-2-gage.mat');

B = load('./allegheny_data/setup/setup-233-gage.mat');
basin = B.basin;
gage = B.gage;
HH = B.HH;

basin.distmat = aa.basin.distmat;
load('./allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('./allegheny_data/setup/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

clearvars aa B

%% Cut down to a shorter time period

i1=60; i2=179;
tv = tv(i1:i2);
gage = gage(i1:i2,:);
nldas_runoff = nldas_runoff(:,:,i1:i2);
true_runoff = true_runoff(:,i1:i2); 
[nt,m] = size(gage);
nldas.runoff  = nldas.runoff(:, :, i1:i2);
tmpa.runoff  = tmpa.runoff(:, :, i1:i2);

%% Assemble runoff prior (Y20)

figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;

tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

%% Assemble runoff prior (Ensemble)

addgauss=1;
multlog=0;

M = 500; % ensemble size

totmean = mean(nldas.runoff_mean);
runoff_init = ones(nt, n)*totmean;

prior_runoff_ens = zeros(n,nt,M);
if addgauss
    
%     a = load('./allegheny_data/errors/m0a1L40T5_AG/Y20_prior_AG.mat');
%     a = load('./allegheny_data/errors/m0a1L40T5_AG/runoff_errors_for_variogram.mat');
%     a = load('./allegheny_data/errors/m0a1L40T5_prop_prior/ens_AG_m0a1L40T5_tmpa_ensemble.mat');
%     a = load('./allegheny_data/errors/m0a1L40T5_prop_prior/ens_AG_m0a1L0T0_tmpa_ensemble.mat');
%       a = load('./allegheny_data/errors/m0a1L40T5_prop_prior/ens_AG_m0a1L40T5_uniform_ensemble.mat');
    a = load('./allegheny_data/errors/m0a1L40T5_prop_prior/ens_AG_m0a1L0T0_uniform_ensemble.mat');
    runoff_errors = a.runoff_errors;
    clearvars a
    %mm=1; % ensemble is centered around Y20 prior
    %Y20_runoff = true_runoff + runoff_errors(:,:,mm); 
    for mm=1:M
        prior_runoff_ens(:,:,mm) = tmpa_runoff_prior' + runoff_errors(:,:,mm);
%         prior_runoff_ens(:,:,mm) = true_runoff + runoff_errors(:,:,mm);
        %prior_runoff_ens(:,:,mm) = Y20_runoff + runoff_errors(:,:,mm);
    end
    
%     a = load('./allegheny_data/errors/m0aflatL40T5_AG/Ens_prior_AGflat.mat');
%     prior_runoff_ens = a.Ens_prior;
%     prior_runoff_ens = a.prior_runoff_ens;
    
elseif multlog
    a = load('./allegheny_data/errors/m1a1L40T5_LM/simulated_runoff_errors_2000.mat');
    runoff_errors = reshape(a.alldata, 20*21, nt, M);
    clearvars a
    basinmask = basin.mask; % no nan values
    basinmask(isnan(basin.mask)) = 0;
    basinmask = logical(basinmask);
    runoff_errors = runoff_errors(basinmask, :, :); 
    mm=1;
    Y20_runoff = true_runoff.*runoff_errors(:,:,mm);
    prior_runoff_ens = zeros(n, nt, M);
    for mm=2:M
        % prior_runoff_ens(:,:,mm) = Y20_runoff.*runoff_errors(:,:,mm);
        prior_runoff_ens(:,:,mm) = true_runoff.*runoff_errors(:,:,mm);
    end
end

% Plot prior maps

figure
t=73:77;
for i=1:3
%     prior_runoff_map_ENS = make_map(basin, squeeze(mean(prior_runoff_ens(:,t(i),:),3)));
    prior_runoff_map_ENS = make_map(basin, squeeze(mean(runoff_errors(:,t(i),:),3)));
    subplot(1,3,i)
    plotraster(basin.lonv, basin.latv, prior_runoff_map_ENS, ['Prior ensemble mean (day ' num2str(t(i)) ')'])
%     caxis([0,6])
end
colormap cool

%% Generate discharge "measurements" (Ensemble)

rng(704753262,'twister')

true_discharge = gage;
swot = 0;
percent_error = 0;

m = size(gage,2);

addmeas_error=1;

if addmeas_error
    disp('using additive meas error')

    mu1 = 0; % mean of additive error
    add_err = zeros(nt,m);
    for tt=1:nt
        for mm=1:m
            sigma1 = (percent_error/100)*gage(tt,mm);
            add_err(tt,mm) = mu1 + sigma1*randn(1,1);
        end
    end
    error_corrupted_discharge_meas = gage + add_err;    
    
else
    % ENKF() assumes homoscedastic, uncorrelated lognormal errors here
    % we will use unbiased relative error of 1% (COV=0.01, mu = 1)
    mean1 = 1;
    Cv = (percent_error/100)^2; % variance
    mu1 = log((mean1^2)/sqrt(Cv+mean1^2));
    sigma1 = sqrt(log(Cv/(mean1^2)+1));
    error_corrupted_discharge_meas = true_discharge.*lognrnd(mu1, sigma1, nt, m);
end

freq = 10;
if swot
    allobs = error_corrupted_discharge_meas;
    swot_obs = NaN(size(allobs));
    swot_obs(1:freq:end,:) = allobs(1:freq:end,:);   
    error_corrupted_discharge_meas = swot_obs;
end

% Plot measurements for each gage
figure
lw = 2;
ind = 1;
for gg=[29]
    subplot(1,1,ind)
    plot(tv, true_discharge(:,gg), 'linewidth', lw)
    hold on
    plot(tv, error_corrupted_discharge_meas(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
    ylim([0,800])
end

%% Set aside validation gauges 

rng(704753262, 'twister');
nval = 7;
val_ind = 1 - 1 + randperm(14-1+1,nval); % setting aside nval gauges for validation
cal_ind = find(~ismember(1:14, val_ind));

% Plot map showing gauge locations

% two gauge case
cal_ind=[29, 146]; % 29, 146
val_ind = [];

% cal_ind = 2;
% val_ind = [1,3:14];

% cal_ind=1:14;
% val_ind = []

ms = 20;
res = 1/8;
figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Allegheny basin')
hold on
plot(basin.gage_lon(cal_ind)-res/2, basin.gage_lat(cal_ind)-res/2, 'r.', 'markersize', ms)
plot(basin.gage_lon(val_ind)-res/2, basin.gage_lat(val_ind)-res/2, 'green.', 'markersize', ms)
legend('', 'Calibration gauges','Validation gauges')

orig_basin = basin;
basin.gage_area = basin.gage_area(cal_ind);
basin.mask_gage = basin.mask_gage(:,:,cal_ind);
basin.flow_length_gage = basin.flow_length_gage(:,:,cal_ind);
basin.gage_i = basin.gage_i(cal_ind);
basin.gage_j = basin.gage_j(cal_ind);
basin.gage_lon = basin.gage_lon(cal_ind);
basin.gage_lat = basin.gage_lat(cal_ind);

%% Do ISR (Ensemble)

%save('~/Downloads/m0s2L40T5_ens_inputs.mat')

% For some reason I need to load this instance of "basin"...
% Probably due to flipud tbh
C = load('./allegheny_data/setup/setup-2-gage.mat');
C.basin.mask = flipud(C.basin.mask);

% Trim down to size for test purposes
% tmax = 50;
% Mmax = 10; 
% prior_runoff_ens = prior_runoff_ens(:,1:tmax,1:Mmax);
% error_corrupted_discharge_meas = error_corrupted_discharge_meas(1:tmax,:);
% % A.basin
% tv = tv(1:tmax);

runoff_prior = ones(n, nt, M);
Cv = (percent_error/100)^2;
s = 2*(k+1); % window length is sensitive is T>window length
C.basin.flow_vel = flow_vel;
C.basin.timestep = 3600*24;
orig_basin = basin;
basin.flow_vel=flow_vel;
basin.timestep=3600*24;

figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask')

opt.gauss = 1; % use additive Gaussian errors (vs. lognormal mult errors)
opt.pbs = 0; % option for particle batch smoother (not working)
opt.multerrors = 0; % using multiplicative (log) update (i.e. discharge errors are mult)
% without localization
opt.loc = 0;
[post_runoff, ~, ~, p] = ISR_EnKF(prior_runoff_ens, ...
    HH(cal_ind,:,:), error_corrupted_discharge_meas(:,cal_ind), s, basin, Cv, opt);

truth = struct();
truth.total_runoff = true_runoff;
truth.true_runoff = true_runoff;
basin.true_runoff = truth.true_runoff;
mean_post_runoff = mean(post_runoff,3);

%%
save('./allegheny_data/results/E5-6/uniform_prior_Cv0/E5_L0T0_both_unif.mat', ...
    'prior_runoff_ens','post_runoff','mean_post_runoff', ...
    'truth', 'error_corrupted_discharge_meas','basin', ...
    'cal_ind', 'val_ind');


% E5_L0T0_outlet
% E5_L0T0_upstream
% E5_L0T0_both
% E6_L40T5_outlet
% E6_L40T5_upstream
% E6_L40T5_both

%% Figure 4.3.1/2

% Compare subbasin cases
gi = (k+1):nt-(k+1);

% Uniform runoff prior, no error
% E6.ds = load('./allegheny_data/results/E5-6/E6_L40T5_outlet_unif.mat');
% E6.us = load('./allegheny_data/results/E5-6/E6_L40T5_upstream_unif.mat');
% E6.both = load('./allegheny_data/results/E5-6/E6_L40T5_both_unif.mat');
% E6.ds = load('./allegheny_data/results/E5-6/uniform_prior_Cv0/E5_L0T0_outlet_unif.mat');
% E6.us = load('./allegheny_data/results/E5-6/uniform_prior_Cv0/E5_L0T0_upstream_unif.mat');
% E6.both = load('./allegheny_data/results/E5-6/uniform_prior_Cv0/E5_L0T0_both_unif.mat');

% TMPA runoff prior
% E6.ds = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E6_L40T5_outlet.mat');
% E6.us = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E6_L40T5_upstream.mat');
% E6.both = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E6_L40T5_both.mat');
E5.ds = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E5_L0T0_outlet.mat');
E5.us = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E5_L0T0_upstream.mat');
E5.both = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E5_L0T0_both.mat');
E6.ds = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E6_L40T5_outlet.mat');
E6.us = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E6_L40T5_upstream.mat');
E6.both = load('./allegheny_data/results/E5-6/tmpa_prior_Cv5/E6_L40T5_both.mat');

% Plot change in NSE from prior to posterior for each case

[~, ~, ~, E5.ds.nsemap] = plot_gofmaps(basin, E5.ds.mean_post_runoff', true_runoff, gi);
[~, ~, ~, E5.us.nsemap] = plot_gofmaps(basin, E5.us.mean_post_runoff', true_runoff, gi);
[~, ~, ~, E5.both.nsemap] = plot_gofmaps(basin, E5.both.mean_post_runoff', true_runoff, gi);
[~, ~, ~, E6.ds.nsemap] = plot_gofmaps(basin, E6.ds.mean_post_runoff', true_runoff, gi);
[~, ~, ~, E6.us.nsemap] = plot_gofmaps(basin, E6.us.mean_post_runoff', true_runoff, gi);
[~, ~, ~, E6.both.nsemap] = plot_gofmaps(basin, E6.both.mean_post_runoff', true_runoff, gi);
[~, ~, ~, nsemap_prior] = plot_gofmaps(basin, tmpa_runoff_prior, true_runoff, gi);

basin.mask = flipud(basin.mask);
E5.ds.postmeanmap = make_map(basin, mean(E5.ds.mean_post_runoff,2));
E5.us.postmeanmap = make_map(basin, mean(E5.us.mean_post_runoff,2));
E5.both.postmeanmap = make_map(basin, mean(E5.both.mean_post_runoff,2));
E6.ds.postmeanmap = make_map(basin, mean(E6.ds.mean_post_runoff,2));
E6.us.postmeanmap = make_map(basin, mean(E6.us.mean_post_runoff,2));
E6.both.postmeanmap = make_map(basin, mean(E6.both.mean_post_runoff,2));
priormeanmap = make_map(basin, mean(tmpa_runoff_prior',2));

figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, flipud(E6.us.postmeanmap - priormeanmap), 'runoff change (us)')
caxis([-.1,0.4])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, flipud(E6.ds.postmeanmap - priormeanmap), 'runoff change (ds)')
caxis([-.1,0.4])
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, flipud(E6.both.postmeanmap - priormeanmap), 'runoff change (both)')
caxis([-.1,0.4])

% Check that runoff matches at the gage
E6.us.postQ = state_model_dumb(E6.us.mean_post_runoff', HH(cal_ind,:,:));
figure
plot(E6.us.postQ(:,2))
hold on
plot(true_discharge(:,146));
legend('Posterior','True')

% plot outline of subdomain (export raster, convert to vector, plot vector)
% R1 = makerefmat(min(basin.lonv), min(basin.latv), 1/8,1/8);
% geotiffwrite('./allegheny_data/results/subdomain1_mask.tif', basin.mask_gage(:,:,1), R1)
% geotiffwrite('./allegheny_data/results/subdomain2_mask.tif', basin.mask_gage(:,:,2), R1)
subdomains_vect = shaperead('./allegheny_data/results/allegheny_subdomains/subdomains.shp');

subdomain1 = E6.both.basin.mask_gage(:,:,1) - E6.both.basin.mask_gage(:,:,2);

%% Get NSE changes for each subdomain

% Method for getting indices of cells in each subdomain
[aa,bb]=(find(subdomain1==1));
subdomain1_index = sub2ind(size(basin.mask), aa,bb);

[aa,bb]=(find(E6.both.basin.mask_gage(:,:,2)==1));
subdomain2_index = sub2ind(size(basin.mask), aa,bb);

% get runoff for grid cells in each subdomain
% Reshape TMPA runoff prior to map

E5.sub1runoff = struct();
E5.sub2runoff = struct();

E5.sub1runoff.prior = get_runoff_for_cells_in_subdomain(tmpa_runoff_prior, subdomain1_index, basin);
E5.sub2runoff.prior = get_runoff_for_cells_in_subdomain(tmpa_runoff_prior, subdomain2_index, basin);
E6.sub1runoff.prior = get_runoff_for_cells_in_subdomain(tmpa_runoff_prior, subdomain1_index, basin);
E6.sub2runoff.prior = get_runoff_for_cells_in_subdomain(tmpa_runoff_prior, subdomain2_index, basin);

E5.sub1runoff.us = get_runoff_for_cells_in_subdomain(E5.us.mean_post_runoff', subdomain1_index, basin);
E5.sub2runoff.us = get_runoff_for_cells_in_subdomain(E5.us.mean_post_runoff', subdomain2_index, basin);
E6.sub1runoff.us = get_runoff_for_cells_in_subdomain(E6.us.mean_post_runoff', subdomain1_index, basin);
E6.sub2runoff.us = get_runoff_for_cells_in_subdomain(E6.us.mean_post_runoff', subdomain2_index, basin);

E5.sub1runoff.ds = get_runoff_for_cells_in_subdomain(E5.ds.mean_post_runoff', subdomain1_index, basin);
E5.sub2runoff.ds = get_runoff_for_cells_in_subdomain(E5.ds.mean_post_runoff', subdomain2_index, basin);
E6.sub1runoff.ds = get_runoff_for_cells_in_subdomain(E6.ds.mean_post_runoff', subdomain1_index, basin);
E6.sub2runoff.ds = get_runoff_for_cells_in_subdomain(E6.ds.mean_post_runoff', subdomain2_index, basin);

E5.sub1runoff.both= get_runoff_for_cells_in_subdomain(E5.both.mean_post_runoff', subdomain1_index, basin);
E5.sub2runoff.both = get_runoff_for_cells_in_subdomain(E5.both.mean_post_runoff', subdomain2_index, basin);
E6.sub1runoff.both= get_runoff_for_cells_in_subdomain(E6.both.mean_post_runoff', subdomain1_index, basin);
E6.sub2runoff.both = get_runoff_for_cells_in_subdomain(E6.both.mean_post_runoff', subdomain2_index, basin);

E5.sub1runoff.true= get_runoff_for_cells_in_subdomain(true_runoff', subdomain1_index, basin);
E5.sub2runoff.true = get_runoff_for_cells_in_subdomain(true_runoff', subdomain2_index, basin);
E6.sub1runoff.true= get_runoff_for_cells_in_subdomain(true_runoff', subdomain1_index, basin);
E6.sub2runoff.true = get_runoff_for_cells_in_subdomain(true_runoff', subdomain2_index, basin);

% Instead of individual cells, let's plot time series for the average of
% each basin

%% Plot subdomains maps

res = 1/8;
figure
subplot(2,1,1)
plotraster(basin.lonv, basin.latv, flipud(subdomain1), 'Downstream subbasin')
hold on 
plot(E6.both.basin.gage_lon(1)-res/2, E6.both.basin.gage_lat(1)-res/2, '.b', 'markersize', ms)
text(E6.both.basin.gage_lon(1), E6.both.basin.gage_lat(1), 'Outlet gage', 'fontsize', fs)
plot(E6.both.basin.gage_lon(2)-res/2, E6.both.basin.gage_lat(2)-res/2, '.b', 'markersize', ms)
text(E6.both.basin.gage_lon(2), E6.both.basin.gage_lat(2)+1/16, 'Upstream gage', 'fontsize', fs)
% plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
% text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
% plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
% text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
xlabel('Lon')
ylabel('Lat')
colorbar off

subplot(2,1,2)
plotraster(basin.lonv, basin.latv, flipud(E6.both.basin.mask_gage(:,:,2)), 'Upstream subbasin')
hold on
plot(E6.both.basin.gage_lon(1)-res/2, E6.both.basin.gage_lat(1)-res/2, '.b', 'markersize', ms)
text(E6.both.basin.gage_lon(1), E6.both.basin.gage_lat(1), 'Outlet gage', 'fontsize', fs)
plot(E6.both.basin.gage_lon(2)-res/2, E6.both.basin.gage_lat(2)-res/2, '.b', 'markersize', ms)
text(E6.both.basin.gage_lon(2), E6.both.basin.gage_lat(2)+1/16, 'Upstream gage', 'fontsize', fs)
% plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
% text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
% plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
% text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
xlabel('Lon')
ylabel('Lat')
colorbar off

%% Figure 4.3.1 (revised)

% Subdomain 1, no corr
figure
subplot(2,2,1)
plot(tv, mean(E5.sub1runoff.prior,1), '-r', 'linewidth', lw)
hold on
grid on
plot(tv, mean(E5.sub1runoff.us,1), '--b', 'linewidth', lw)
plot(tv, mean(E5.sub1runoff.ds,1), '--g', 'linewidth', lw)
plot(tv, mean(E5.sub1runoff.both,1), '--cyan', 'linewidth', lw)
plot(tv, mean(E5.sub1runoff.true,1), '-k', 'linewidth', lw)
title('Downstream subbasin mean runoff (L0T0)')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)
axis tight
xlim([tv(20),tv(60)])

% Subdomain 2, no corr
subplot(2,2,2)
plot(tv, mean(E5.sub2runoff.prior,1), '-r', 'linewidth', lw)
hold on
grid on
plot(tv, mean(E5.sub2runoff.us,1), '--b', 'linewidth', lw)
plot(tv, mean(E5.sub2runoff.ds,1), '--g', 'linewidth', lw)
plot(tv, mean(E5.sub2runoff.both,1), '--cyan', 'linewidth', lw)
plot(tv, mean(E5.sub2runoff.true,1), '-k', 'linewidth', lw)
title('Upstream subbasin mean runoff (L0T0)')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)
axis tight
xlim([tv(20),tv(60)])

% Subdomain 1, corr
subplot(2,2,3)
plot(tv, mean(E6.sub1runoff.prior,1), '-r', 'linewidth', lw)
hold on
grid on
plot(tv, mean(E6.sub1runoff.us,1), '--b', 'linewidth', lw)
plot(tv, mean(E6.sub1runoff.ds,1), '--g', 'linewidth', lw)
plot(tv, mean(E6.sub1runoff.both,1), '--cyan', 'linewidth', lw)
plot(tv, mean(E6.sub1runoff.true,1), '-k', 'linewidth', lw)
title('Downstream subbasin mean runoff (L40T5)')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)
axis tight
xlim([tv(20),tv(60)])

% Subdomain 2, corr
subplot(2,2,4)
plot(tv, mean(E6.sub2runoff.prior,1), '-r', 'linewidth', lw)
hold on
grid on
plot(tv, mean(E6.sub2runoff.us,1), '--b', 'linewidth', lw)
plot(tv, mean(E6.sub2runoff.ds,1), '--g', 'linewidth', lw)
plot(tv, mean(E6.sub2runoff.both,1), '--cyan', 'linewidth', lw)
plot(tv, mean(E6.sub2runoff.true,1), '-k', 'linewidth', lw)
title('Upstream subbasin mean runoff (L40T5)')
ylabel('Runoff (mm/day)')
xlabel('Date')
set(gca, 'fontsize', fs)
axis tight
xlim([tv(20),tv(60)])
legend('Prior','Posterior (ups. gage only) ','Posterior (downstr gage only)',...
    'Posterior (both gages)','Truth', 'fontsize', 12)

% % Plot subdomains maps
% 
% res = 1/8;
% 
% subplot(2,3,3)
% plotraster(basin.lonv, basin.latv, flipud(subdomain1), 'Downstream subbasin')
% hold on 
% plot(E6.both.basin.gage_lon(1)-res/2, E6.both.basin.gage_lat(1)-res/2, '.b', 'markersize', ms)
% text(E6.both.basin.gage_lon(1), E6.both.basin.gage_lat(1), 'Outlet gage', 'fontsize', fs)
% plot(E6.both.basin.gage_lon(2)-res/2, E6.both.basin.gage_lat(2)-res/2, '.b', 'markersize', ms)
% text(E6.both.basin.gage_lon(2), E6.both.basin.gage_lat(2)+1/16, 'Upstream gage', 'fontsize', fs)
% % plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
% % text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
% % plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
% % text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
% xlabel('Lon')
% ylabel('Lat')
% colorbar off
% 
% subplot(2,3,6)
% plotraster(basin.lonv, basin.latv, flipud(E6.both.basin.mask_gage(:,:,2)), 'Upstream subbasin')
% hold on
% plot(E6.both.basin.gage_lon(1)-res/2, E6.both.basin.gage_lat(1)-res/2, '.b', 'markersize', ms)
% text(E6.both.basin.gage_lon(1), E6.both.basin.gage_lat(1), 'Outlet gage', 'fontsize', fs)
% plot(E6.both.basin.gage_lon(2)-res/2, E6.both.basin.gage_lat(2)-res/2, '.b', 'markersize', ms)
% text(E6.both.basin.gage_lon(2), E6.both.basin.gage_lat(2)+1/16, 'Upstream gage', 'fontsize', fs)
% % plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
% % text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
% % plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
% % text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
% xlabel('Lon')
% ylabel('Lat')
% colorbar off

myNSE(mean(E5.sub1runoff.true,1), mean(E5.sub1runoff.prior,1))
myNSE(mean(E5.sub1runoff.true,1), mean(E5.sub1runoff.us,1))
myNSE(mean(E5.sub1runoff.true,1), mean(E5.sub1runoff.ds,1))
myNSE(mean(E5.sub1runoff.true,1), mean(E5.sub1runoff.both,1))

myNSE(mean(E5.sub2runoff.true,1), mean(E5.sub2runoff.prior,1))
myNSE(mean(E5.sub2runoff.true,1), mean(E5.sub2runoff.us,1))
myNSE(mean(E5.sub2runoff.true,1), mean(E5.sub2runoff.ds,1))
myNSE(mean(E5.sub2runoff.true,1), mean(E5.sub2runoff.both,1))

myNSE(mean(E6.sub1runoff.true,1), mean(E6.sub1runoff.prior,1))
myNSE(mean(E6.sub1runoff.true,1), mean(E6.sub1runoff.us,1))
myNSE(mean(E6.sub1runoff.true,1), mean(E6.sub1runoff.ds,1))
myNSE(mean(E6.sub1runoff.true,1), mean(E6.sub1runoff.both,1))

myNSE(mean(E6.sub2runoff.true,1), mean(E6.sub2runoff.prior,1))
myNSE(mean(E6.sub2runoff.true,1), mean(E6.sub2runoff.us,1))
myNSE(mean(E6.sub2runoff.true,1), mean(E6.sub2runoff.ds,1))
myNSE(mean(E6.sub2runoff.true,1), mean(E6.sub2runoff.both,1))

%%
% Calculate average NSE for cells in each subdomain

[E5.sub1runoff, E5.sub2runoff] = calc_med_runoff_nse_for_localiz_fig(E5.sub1runoff, E5.sub2runoff, gi);
[E6.sub1runoff, E6.sub2runoff] = calc_med_runoff_nse_for_localiz_fig(E6.sub1runoff, E6.sub2runoff, gi);

%% Figure 4.3.1

cmin = -.0;
cmax = 1;
fmtstr = 'blue-';
col1 = 'white';
col2 = 'black';
lw2 = 4;

figure
% plotraster(basin.lonv, basin.latv, flipud(nsemap_us), 'Upstream')
plotraster(basin.lonv, basin.latv, (nsemap_prior), 'Prior')
% plotraster(basin.lonv, basin.latv, (E6.us.nsemap - nsemap_prior), 'Upstream')
% plotraster(basin.lonv, basin.latv, flipud(E6.us.postmeanmap - priormeanmap), 'runoff change (us)')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
% text(-79.75,41.25,'\DeltaNSE = 0', 'color', col2, 'fontsize', fs)
% text(-79.5+2/4,42,'\DeltaNSE = 1.1', 'color', col2, 'fontsize', fs)
plot(E6.both.basin.gage_lon - res/2, E6.both.basin.gage_lat - res/2, '.k', 'markersize', ms)
caxis([cmin,cmax])

figure
subplot(1,3,1)
% plotraster(basin.lonv, basin.latv, flipud(nsemap_us), 'Upstream')
plotraster(basin.lonv, basin.latv, (E6.us.nsemap), 'Upstream')
% plotraster(basin.lonv, basin.latv, (E6.us.nsemap - nsemap_prior), 'Upstream')
% plotraster(basin.lonv, basin.latv, flipud(E6.us.postmeanmap - priormeanmap), 'runoff change (us)')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
% text(-79.75,41.25,'\DeltaNSE = 0', 'color', col2, 'fontsize', fs)
% text(-79.5+2/4,42,'\DeltaNSE = 1.1', 'color', col2, 'fontsize', fs)
plot(E6.both.basin.gage_lon - res/2, E6.both.basin.gage_lat - res/2, '.k', 'markersize', ms)
caxis([cmin,cmax])
title('Upstream gage only')

subplot(1,3,2)
% plotraster(basin.lonv, basin.latv, flipud(nsemap_ds), 'Upstream')
plotraster(basin.lonv, basin.latv, (E6.ds.nsemap), 'Downstream')
% plotraster(basin.lonv, basin.latv, (E6.ds.nsemap - nsemap_prior), 'Downstream')
% plotraster(basin.lonv, basin.latv, flipud(E6.ds.postmeanmap - priormeanmap), 'runoff change (ds)')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
% % text(-79.75,41.25,'\DeltaNSE = 0.87', 'color', col2, 'fontsize', fs)
% text(-79.5+2/4,42,'\DeltaNSE = 0.94', 'color', col2, 'fontsize', fs)
plot(E6.both.basin.gage_lon - res/2, E6.both.basin.gage_lat - res/2, '.k', 'markersize', ms)
caxis([cmin,cmax])
title('Downstream gage only')

subplot(1,3,3)
% plotraster(basin.lonv, basin.latv, flipud(nsemap_both), 'Upstream')
plotraster(basin.lonv, basin.latv, (E6.both.nsemap), 'Both')
% plotraster(basin.lonv, basin.latv, flipud(E6.both.postmeanmap - priormeanmap), 'runoff change (both)')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
% text(-79.75,41.25,'\DeltaNSE = 0.89', 'color', col2, 'fontsize', fs)
% text(-79.5+2/4,42,'\DeltaNSE = 1.1', 'color', col2, 'fontsize', fs)
plot(E6.both.basin.gage_lon - res/2, E6.both.basin.gage_lat - res/2, '.k', 'markersize', ms)
caxis([cmin,cmax])
title('Both gages')
% colormap(colorMap)
% colormap(flipud(cool))
% colormap(bluewhitered)

%% Show differences between each case

figure
subplot(1,2,1)
% plotraster(basin.lonv, basin.latv, flipud(nsemap_both), 'Upstream')
plotraster(basin.lonv, basin.latv, (E6.both.nsemap - E6.ds.nsemap), 'both - ds')
% plotraster(basin.lonv, basin.latv, flipud(E6.both.postmeanmap - E6.ds.postmeanmap), '\Delta runoff both-ds')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
% text(-79.75,41.25,'\DeltaNSE = 0.89', 'color', col2, 'fontsize', fs)
% text(-79.5+2/4,42,'\DeltaNSE = 1.1', 'color', col2, 'fontsize', fs)
plot(E6.both.basin.gage_lon - res/2, E6.both.basin.gage_lat - res/2, '.k', 'markersize', ms)
caxis([cmin,cmax])

subplot(1,2,2)
% plotraster(basin.lonv, basin.latv, flipud(nsemap_both), 'Upstream')
plotraster(basin.lonv, basin.latv, (E6.both.nsemap - E6.us.nsemap), 'both - us')
% plotraster(basin.lonv, basin.latv, flipud(E6.both.postmeanmap - E6.us.postmeanmap), '\Delta runoff both-us')
hold on
plot(subdomains_vect(1).X, subdomains_vect(1).Y, fmtstr, 'linewidth', lw2)
plot(subdomains_vect(2).X, subdomains_vect(2).Y, fmtstr, 'linewidth', lw2)
% text(-79.75,41.25,'\DeltaNSE = 0.89', 'color', col2, 'fontsize', fs)
% text(-79.5+2/4,42,'\DeltaNSE = 1.1', 'color', col2, 'fontsize', fs)
plot(E6.both.basin.gage_lon - res/2, E6.both.basin.gage_lat - res/2, '.k', 'markersize', ms)
caxis([cmin,cmax])

colormap(bluewhitered)

subdomain2 = flipud(E6.both.basin.mask_gage(:,:,2));
sub1mask = subdomain1==1;
sub2mask = subdomain2==1;

median(E6.both.nsemap(flipud(sub1mask)) - E6.us.nsemap(flipud(sub1mask))) % marginal benefit to sub1 for adding ds gage
median(E6.both.nsemap(sub2mask) - E6.us.nsemap(sub2mask)) % marginal benefit to sub2 for adding ds gage

median(E6.both.nsemap(flipud(sub1mask)) - E6.ds.nsemap(flipud(sub1mask))) % marginal benefit to sub1 for adding us gage
median(E6.both.nsemap(sub2mask) - E6.ds.nsemap(sub2mask)) % marginal benefit to sub2 for adding us gage

%% Figure 4.3.3 EnISR outputs with localization



