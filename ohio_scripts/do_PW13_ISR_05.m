%% PW13 ISR
%
% 6/5/2023 JRS

clear, clc, close all
cd /hdd/ISR/inverse_streamflow_routing
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

load('./ohio_data/ohio_tmpa_nanfilled.mat')
load('./ohio_data/ohio_nldas_nanfilled.mat')
load('./ohio_data/swot_like_measurements_100m_no_error_revised.mat')

% add distmat
A = load('./ohio_data/setup-1-gage.mat');
basin.distmat = A.basin.distmat;

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

% load('./basin_setup/basin.mat')

%% If desired, crop to a shorter temporal domain

% tmax = 50;
% tv = tv(1:tmax);
% true_discharge = true_discharge(1:tmax,:);
% tmpa.runoff = tmpa.runoff(:,:,1:tmax);
% nldas.runoff = nldas.runoff(:,:,1:tmax);

[nt,m] = size(true_discharge);

%% Plot tmpa prior vs. nldas (true) runoff (Y20 Figure 6)

% TMPA prior, NLDAS truth

figure
t=75:77;
cmax = 2.5;
for i=1:3
    subplot(2,3,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
    subplot(2,3,3+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
end
colormap cool

%% Assemble runoff prior (PW13)

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

%% Generate discharge "measurements" (PW13)

% gage = true_discharge; % should corrupt with error

gage_all = true_discharge;
gage = true_discharge; % should corrupt with error
% gage = true_discharge_w_swot_sampling; % should corrupt with error

% Y20: additive Gaussian error with mu, sigma

mu1 = 0; % mean of additive error
sigma1 = 0.30*gage; % stddev of additive error (30% of truth)
add_err = mu1 + sigma1.*randn(nt,m);
gage_w_error = gage + add_err;

% Plot error for several gages
figure
lw = 2;
ind = 1;
for gg=[1,5,10]
    subplot(1,3,ind)
    plot(tv, gage_all(:,gg), 'linewidth', lw)
    hold on
    plot(tv, gage_w_error(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end

err = true_discharge(:) - gage_w_error(:);
nanmean(err)
nanstd(err)

figure
histogram(err)
% why does this histogram not look gaussian?
% because each outcome is from a different distribution, with a different
% standard deviation (sigma varies for each observation)
% If we repeated this many times, then there would be Gaussian dists

% y = rand(1000,50);
% yobs = y + randn(1000, 50);
% figure
% histogram(yobs-y);

%% Plot discharge "measurements" (PW13)

ms = 30;
fs = 16;
figure
plot(tv, true_discharge(:,1), 'linewidth', lw)
hold on 
plot(tv, gage_w_error(:,1), '.', 'markersize', ms)
legend('Truth','Measurements')
xlabel('time')
ylabel('Q (mm/day)')
grid on
title('Discharge at outlet')
set(gca, 'fontsize', fs)

% true_discharge_w_swot_sampling

%% Do ISR (PW13)

% runoff_init = ones(nt,n);

s = k+1;
% s = 2*(k+1)-1; % for window of length 32 days
cov = 1; % coefficient of variation
% R = (0.15^2); % meas error covariance
R = 0;

totmean = mean(nldas.runoff_mean);
totsd = std(nldas.runoff_mean);
runoff_init = ones(nt, n)*totmean;

% takes about 1.5 hr with these parameters on my office PC
tic
% [post_runoff_PW13] = ISR_PW13(tmpa_runoff_prior, HH, gage, s, 'proportional', cov, R);
[post_runoff_PW13] = ISR_PW13(runoff_init, HH, gage, s, 'const_diag', totsd, R);
toc

save('./ohio_data/ISR_results_PW13_m240_daily_no_err_null.mat', 'post_runoff_PW13', 's', 'cov', 'R', 'gage')

%% Plot overview of results

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;

gi = (k+1):nt-(k+1);
plt_ISR_results_overview(basin, tmpa_runoff_prior', post_runoff_PW13', truth, tv, gi)

