%% Comparing Y20 ISR and ensemble ISR
%
% 6/5/2023 JRS

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

load('./ohio_data/ohio_tmpa_nanfilled.mat')
load('./ohio_data/ohio_nldas_nanfilled.mat')
load('./ohio_data/swot_like_measurements_100m_no_error.mat')

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
t=43:47;
for i=1:5
    subplot(2,5,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
    caxis([0,6])
    subplot(2,5,5+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
    caxis([0,6])
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

gage = true_discharge; % should corrupt with error

% gage = true_discharge_w_swot_sampling; % should corrupt with error

% Y20: additive Gaussian error with mu, sigma

mu1 = 0; % mean of additive error
% sigma1 = 0.3*gage; % stddev of additive error (30% of truth)
for tt=1:nt
    for mm=1:m
        sigma1 = 0.15*gage(tt,mm);
        add_err(tt,mm) = mu1 + sigma1*randn(1,1);
    end
    disp(tt)
end
gage_w_error = gage + add_err;
% add_err = mu1 + sigma1.*randn(nt,m);

% Plot error for several gages
figure
lw = 2;
ind = 1;
for gg=[1,5,10]
    subplot(1,3,ind)
    plot(tv, gage(:,gg), 'linewidth', lw)
    hold on
    plot(tv, gage_w_error(:,gg), 'red.', 'markersize', 20)
    xlabel('Time')
    legend('Actual discharge','Measured discharge')
    ylabel('Discharge (mmd/day)')
    ind = ind + 1;
end
% ensemble ISR: multiplicative LN error with mean, stddev

err = true_discharge(:) - gage_w_error(:);
nanmean(err)
nanstd(err)

figure
histogram(err)

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

s = 2*(k+1);
alpha1 = 1;
R = (0.15^2);

% takes about 1.5 hours with these parameters on my PC
tic
[post_runoff_PW13, Klast] = ISR_PW13(tmpa_runoff_prior, HH, gage, s, 'proportional', alpha1, R);
toc
% COV = 1 and COV = 10 do not produce dramatically different estimates

%% Do ISR (Ensemble)

% For some reason I need to load this instance of "basin"...
% Probably due to flipud tbh
C = load('./basin_setup/swot_like_measurements_100m_no_error.mat');

% Trim down to size for test purposes
% tmax = 50;
% Mmax = 10; 
% prior_runoff_ens = prior_runoff_ens(:,1:tmax,1:Mmax);
% error_corrupted_discharge_meas = error_corrupted_discharge_meas(1:tmax,:);
% % A.basin
% tv = tv(1:tmax);

% Could the issue be with make_map/flipud? Could that be why we are getting
% NaN/Inf?
opt = struct();
runoff_prior = ones(n, nt, M);
Cv = (0.01)^2;
s = k+1;
C.basin.flow_vel = flow_vel;
C.basin.timestep = 3600*24;
tic
[post_runoff, w1, w2, maxv, maxv_sub] = ISR_EnKF_domain_wrapper_d2(prior_runoff_ens, ...
    HH, error_corrupted_discharge_meas, s, C.basin, Cv, opt, tv);
toc

% Save results
save('./ohio_basin_ensemble_ISR_results_M200_Cv0.01.mat', 'post_runoff', '-v7.3')

% post_runoff_orig = post_runoff;
% post_runoff(post_runoff>1e5) = 0;
% post_runoff(isnan(post_runoff)) = 0;

% post_runoff_orig = post_runoff;
% imag_ind = find(abs(imag(post_runoff))>0);
% post_runoff(imag_ind) = NaN;

%%

ENS = struct(); % structure to store results
ENS.mean_prior_runoff = mean(prior_runoff_ens,3);       
ENS.sd_prior_runoff = std(prior_runoff_ens, [], 3);
ENS.mean_posterior_runoff = mean(post_runoff, 3);
ENS.sd_posterior_runoff = std(post_runoff, [], 3);

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;
gi = (k+1):nt-(k+1);
plot_ensemble_runoff_maps(basin, ENS.mean_prior_runoff, ENS.mean_posterior_runoff, truth, tv, gi)

%% Load ensemble ISR results (use ohio_domain_ISR_script)

dom = load('./domain_ISR_results_M5.mat');

%% Load previous PW13/Y20 ISR results

Y20 = load('./ohio_results_PW13_Y20.mat');

%% Compare prior, posterior, truth (basin mean)

tv = datetime(2009,1,1):datetime(2009,12,31);
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