%% Comparing Y20 ISR and ensemble ISR
%
% 7/21/2023 JRS
% Implementing domain Y20 ISR for Allegheny basin

clear, clc, close all
addpath(genpath('./src/'))

A=load('./allegheny_data/setup/setup-2-gage.mat');
load('./allegheny_data/setup/setup-swot-gage.mat');
basin.distmat = A.basin.distmat; clearvars A
load('./allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('./allegheny_data/setup/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);
true_discharge = gage;
[nt,m] = size(true_discharge);

% basin.mask = flipud(basin.mask);

%% Cut down to a shorter time period

i1=60; i2=179;
% i1=60; i2=150;
tv = tv(i1:i2);
gage = gage(i1:i2,:);
nldas_runoff = nldas_runoff(:,:,i1:i2);
tmpa_runoff = tmpa_runoff(:,:,i1:i2);
true_runoff = true_runoff(:,i1:i2); 
[nt,m] = size(gage);
tmpa.runoff = tmpa.runoff(:, :, i1:i2);
nldas.runoff  = nldas.runoff(:, :, i1:i2);
true_discharge = true_discharge(i1:i2,:);
nt=length(tv);

%% Plot tmpa prior vs. nldas (true) runoff (Y20 Figure 6)

% TMPA prior, NLDAS truth

figure
t=43:47;
for i=1:2
    subplot(2,2,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
    caxis([0,6])
    subplot(2,2,2+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
    caxis([0,6])
end
colormap cool

%% Assemble runoff prior (Y20)

% basin.mask = flipud(basin.mask);
figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;
tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

%% Generate a runoff prior with known error characteristics

% This method is RAM-limited bc errors are generated in one batch
M=1000; % number of replicates

% Parameters for generating runoff error covariance
opt.L = 40; % km
opt.T = 5; % days 
opt.RAM = 36;
opt.rho_thres = 0;
opt.progress=1;
opt.alpha = 2; % proportionality constant
opt.SC_savename = './sc.mat';
opt.TC_savename = './tc.mat';
opt.rho_savename = './rho.mat';

% T = 5, 25, 1e6 discharge estimates comparison...

% Generate normally-distributed errors
mu1 = zeros(n*nt,1); % mean of runoff errors
x4cov = fliparoo(true_runoff(1:n,1:nt)', nt-1); % ok
sigma1 = calc_yang_covariance2(basin, opt, n, nt-1, x4cov); % covariance of runoff errors
tmp = mvnrnd(mu1, sigma1, M)'; % runoff errors
runoff_errors = reshape(tmp, n, nt, M); % ok

save('./allegheny_data/errors/m0s2L40T5_norm_covmat.mat','sigma1')

% Generate prior
prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
    prior_runoff_ens(:,:,mm) = true_runoff + runoff_errors(:,:,mm);
%     prior_runoff_ens(:,:,mm) = tmpa_runoff_prior'.*lognrnd(mu1, sigma1, n, nt);
end
disp('Generated prior ensemble')

%% Generate discharge "measurements" (Y20)

rng(704753262,'twister')

% gage = true_discharge; % should corrupt with error
% gage = true_discharge_w_swot_sampling; % should corrupt with error

% Y20: additive Gaussian error with mu, sigma
swot = 1;
freq = 10;
percent_error = 15;

mu1 = 0; % mean of additive error
% sigma1 = 0.3*gage; % stddev of additive error (30% of truth)
add_err = zeros(nt, m);
for tt=1:nt
    for mm=1:m
        sigma1 = (percent_error/100)*gage(tt,mm);
        add_err(tt,mm) = mu1 + sigma1*randn(1,1);
    end
    disp(tt)
end
gage_w_error = gage + add_err;

if swot  
   tmp = gage_w_error;
   gage_w_error = NaN(size(tmp));
   gage_w_error(1:freq:end,:) = tmp(1:freq:end,:);
end

% add_err = mu1 + sigma1.*randn(nt,m);

% Plot error for several gages
figure
lw = 2;
ind = 1;
for gg=[2,14]
    subplot(1,2,ind)
    plot(tv, true_discharge(:,gg), 'linewidth', lw)
    hold on
    plot(tv, gage_w_error(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end
% ensemble ISR: multiplicative LN error with mean, stddev

err = true_discharge(:) - gage_w_error(:);
mean(err)
std(err)

figure
histogram(err)


%% Set aside validation gauges 

rng(704753262, 'twister');
nval = 7; 
val_ind = 1 - 1 + randperm(14-1+1,nval); % setting aside nval gauges for validation
cal_ind = find(~ismember(1:14, val_ind));

cal_ind = 1:14;
val_ind = [];

% cal_ind = [235,240];
% val_ind = 1:240;
% val_ind(cal_ind) = [];
% Plot map showing gauge locations

ms = 10;
res = 1/8;
figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'Ohio basin')
hold on
plot(basin.gage_lon(cal_ind)-res/2, basin.gage_lat(cal_ind)-res/2, 'r.', 'markersize', ms)
plot(basin.gage_lon(val_ind)-res/2, basin.gage_lat(val_ind)-res/2, 'green.', 'markersize', ms)
legend('', 'Calibration gauges','Validation gauges')

%% Plot discharge "measurements" (Y20)

gg=2;
ms = 30;
fs = 16;
figure
plot(tv, true_discharge(:,gg), 'linewidth', lw)
hold on 
plot(tv, gage_w_error(:,gg), '.', 'markersize', ms)
legend('Truth','Measurements')
xlabel('time')
ylabel('Q (mm/day)')
grid on
title('Discharge at outlet')
set(gca, 'fontsize', fs)

% true_discharge_w_swot_sampling

%% Do ISR (Y20)

% 173 seconds per window to calculate Cx -> about 16 hours for 334 windows
% Approx. 180 seconds per window to calculate K using pinv or \ (use pinv)
% It ends up taking about 6 minutes per window to run Y20 ISR for Ohio

% basin.mask = flipud(basin.mask);
% runoff_init = ones(nt,n);

s = 2*(k+1);
optionsfile = './allegheny_data/setup/options_alleg.txt';
% rho_thres = exp(-2);
rho_thres = 0;
L = 40; % km
T = 5; % days
alpha1 = 1;

% total_mean = mean(true_runoff(:));
% null_runoff_prior = total_mean*ones(nt,n);

tic
% [post_runoff_Y20_us_gage, ~] = ISR_Y20(tmpa_runoff_prior, HH(2,:,:), gage_w_error(:,2), s, basin, optionsfile, L, T, rho_thres);
% [post_runoff_Y20_both_gage, ~] = ISR_Y20(tmpa_runoff_prior, HH, gage_w_error, s, basin, optionsfile, L, T, rho_thres);

% TMPA prior
[post_runoff_Y20, post_stddev, prior_stddev, ~] = ISR_Y20(tmpa_runoff_prior, HH(cal_ind,:,:), ...
    gage_w_error(:,cal_ind), s, basin, optionsfile, L, T, rho_thres, alpha1);

toc

% Prior with known error characteristics
% post_runoff_Y20 = zeros(nt,n,M);
% for mm=1:M % loop over replicates
%     post_runoff_Y20(:,:,mm) = ISR_Y20(prior_runoff_ens(:,:,mm)', HH(cal_ind,:,:), ...
%         gage_w_error(:,cal_ind), s, basin, optionsfile, L, T, rho_thres, alpha1);
% end

% Save results
% save('./allegheny_data/results/tmpa_prior/Y20_m0a1L40T5.mat', ...
%     'tmpa_runoff_prior','post_runoff_Y20', 'gage_w_error','-v7.3')

save('./allegheny_data/results/tmpa_prior/Y20_m0a1L40T5.mat', '-v7.3')

%% Look at discharge results

runoff_prior = tmpa_runoff_prior;
% runoff_prior = prior_runoff_ens(:,:,mm)';
post_runoff = post_runoff_Y20;

clc
gi=1:nt;
post_discharge = state_model_dumb(post_runoff, HH);
prior_discharge = state_model_dumb(runoff_prior, HH);

gagei = cal_ind(2);

figure
kge_prior = plot_discharge_ts(prior_discharge(:,gagei), true_discharge(:,gagei))

figure
plot(prior_discharge(gi,gagei), 'b-*', 'linewidth', lw)
hold on
kge_post = plot_discharge_ts(post_discharge(gi,gagei), true_discharge(gi,gagei))
% title('Discharge (M=200, 10 days, no err)')
plot(gage_w_error(gi,gagei), 'k.', 'markersize', 20)
% xlim([100,160])
ylim([0,1000])
legend('Prior','Truth','Posterior','Observations')

% pct neg runoff
100*sum(post_runoff(:)<0)/numel(post_runoff) % percent <0

% Calculate median runoff NSE
priornse = zeros(n,1);
postnse = zeros(n,1);
for kk=1:n
    priornse(kk) = myNSE(true_runoff(kk,gi)', tmpa_runoff_prior(gi,kk));
    postnse(kk) = myNSE(true_runoff(kk,gi)', post_runoff(gi,kk));
end
disp(['median NSE: ' num2str(median(priornse))])
disp(['median NSE: ' num2str(median(postnse))])

%% Plot overview of results

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;

gi = (k+1):nt-(k+1);
% plt_ISR_results_overview(basin, null_runoff_prior', post_runoff_Y20', truth, tv, gi)
plt_ISR_results_overview(basin, tmpa_runoff_prior', post_runoff_Y20', truth, tv, gi)

sum(post_runoff_Y20(:)<0)
figure,imagesc(post_runoff_Y20<0), title('cells with negative runoff')

%% Show where negative runoff values happen

t_ind = 169-i1+1;

runoff_map_w_neg = make_map(basin, post_runoff_Y20(t_ind,:));
tr_map = make_map(basin, true_runoff(:,t_ind)');
prior_map = make_map(basin, tmpa_runoff_prior(t_ind,:));
tv(t_ind)

min(post_runoff_Y20(t_ind,:))

figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, prior_map, 'Prior runoff (mm/day)')
caxis([-5,10])
subplot(1,3,2)
plotraster(basin.lonv, basin.latv, runoff_map_w_neg, 'Posterior runoff (mm/day)')
caxis([-5,10])
subplot(1,3,3)
plotraster(basin.lonv, basin.latv, tr_map, 'True runoff (mm/day)')
caxis([-5,10])
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
colormap(colorMap)

%% Look at prior and posterior runoff for a grid cell

% kk=20; % an average cell (20)
kk=134; % a good cell (133)
% kk=50; % a bad cell (50)
% 113   114   115   134 % cells with negative runoff

figure
h1 = plot(tmpa_runoff_prior(:,kk), 'blue', 'linewidth', lw);
hold on
plot(tmpa_runoff_prior(:,kk) + 2*prior_stddev(:,kk), 'blue--', 'linewidth', lw)
plot(tmpa_runoff_prior(:,kk) - 2*prior_stddev(:,kk), 'blue--', 'linewidth', lw)
h2 = plot(post_runoff_Y20(:,kk), 'r', 'linewidth', lw);
plot(post_runoff_Y20(:,kk) + 2*post_stddev(:,kk), 'r--', 'linewidth', lw);
plot(post_runoff_Y20(:,kk) - 2*post_stddev(:,kk), 'r--', 'linewidth', lw);
h3 = plot(nldas_runoff_true(:,kk), 'k', 'linewidth', lw);
legend([h1,h2,h3],'Prior mean runoff','Posterior mean runoff', 'True runoff')
% xlim([20,60])
disp('runoff NSE for this cell:')
myNSE(nldas_runoff_true(:,kk), tmpa_runoff_prior(:,kk))
myNSE(nldas_runoff_true(:,kk), post_runoff_Y20(:,kk))

%% Save results

