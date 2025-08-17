%% Domain ISR with log transform and ensemble Kalman update
%
% 7/20/2023 JRS

clear, clc, close all
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
% cd /hdd/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))
% addpath('/Users/jschap/Documents/MATLAB/from_file_exchange/raacampbell-shadedErrorBar-e0b2e47')

fs=18;lw=2;
load('./allegheny_data/setup/setup-swot-gage.mat');
aa=load('./allegheny_data/setup/setup-2-gage.mat');
basin.distmat = aa.basin.distmat;
load('./allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('./allegheny_data/setup/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

%% Cut down to a shorter time period

i1=60; i2=179;
% i1=60; i2=90;
tv = tv(i1:i2);
gage = gage(i1:i2,:);
nldas_runoff = nldas_runoff(:,:,i1:i2);
tmpa_runoff = tmpa_runoff(:,:,i1:i2);
true_runoff = true_runoff(:,i1:i2); 
[nt,m] = size(gage);
tmpa.runoff = tmpa.runoff(:, :, i1:i2);
nldas.runoff  = nldas.runoff(:, :, i1:i2);

%% Plot tmpa prior vs. nldas (true) runoff (Y20 Figure 6)

% TMPA prior, NLDAS truth

figure
t=[20,45];
for i=1:2
    subplot(2,2,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
%     caxis([0,6])
    subplot(2,2,2+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
%     caxis([0,6])
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

%% Assemble runoff prior (Ensemble)

% Load in the cosmos errors and use them to corrupt the true runoff
% OR
% Corrupt prior with uncorrelated errors generated using lognrnd

tic
M = 2000; % ensemble size
% A = load('./allegheny_data/errs/m1s2L40T5_lognorm/simulated_runoff_errors_500.mat');
% A = load('./allegheny_data/m0s2L40T5_norm/simulated_runoff_errors_200.mat');
% A = load('./allegheny_data/alleg_errs_m1cv1_L133_T15.8/simulated_runoff_errors_200.mat');
% A = load('./allegheny_data/alleg_errs_m1cv1L0T0_M1000/simulated_runoff_errors_200.mat');

addgauss_errors = 1;
if addgauss_errors
    A = load('./allegheny_data/errors/m0a1L40T5/Ens_errors_AG.mat');
    runoff_errors = A.runoff_errors;
else
    
    % A = load('./allegheny_data/errors/m1a1L40T5_LM/simulated_runoff_errors_2000.mat');

    m1 = mean(A.alldata(:));
    std(A.alldata(:))

    % mu1 = mean(log(A.alldata(:)));
    % std(log(A.alldata(:)))

    A.alldata = A.alldata(:,:,1:nt,:);

    aa=A.alldata(1,1,:,1);
    figure,plot(aa(:))

    % figure
    % histogram(log(A.alldata(:)))
    % % xlim([0,6])
    % xlabel('log Runoff error (mm/day)')
    % ylabel('Count')

    runoff_errors = reshape(A.alldata, 20*21, nt, M);
    basinmask = basin.mask; % no nan values
    basinmask(isnan(basin.mask)) = 0;
    basinmask = logical(basinmask);
    runoff_errors = runoff_errors(basinmask, :, :);

end

total_mean = mean(true_runoff(:));
null_runoff_prior = total_mean*ones(nt,n);

% Corrupt TMPA runoff with random errors to generate prior
prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
%     prior_runoff_ens(:,:,mm) = null_runoff_prior' + runoff_errors(:,:,mm);
    prior_runoff_ens(:,:,mm) = tmpa_runoff_prior' + runoff_errors(:,:,mm);
%     prior_runoff_ens(:,:,mm) = tmpa_runoff_prior'.*runoff_errors(:,:,mm);
end
disp('Generated prior ensemble')
toc

opt.mu1 = mean(log(runoff_errors(:))); % used for bias correction

% Plot prior maps

figure
t=73:77;
for i=1:3
    prior_runoff_map_ENS = make_map(basin, squeeze(mean(runoff_errors(:,t(i),:),3)));
%     prior_runoff_map_ENS = make_map(basin, squeeze(mean(prior_runoff_ens(:,t(i),:),3)));
    subplot(1,3,i)
    plotraster(basin.lonv, basin.latv, prior_runoff_map_ENS, ['Prior ensemble mean (day ' num2str(t(i)) ')'])
%     caxis([0,6])
end
colormap cool

%% Plot time series for two particular grid cells

fs=18;lw=2;

basinmask = flipud(basin.mask);
basinmask(isnan(basinmask)) = 0;
basinmask = logical(basinmask);
% % lon = basin.grid_lon(basinmask);
% lat = basin.grid_lat(basinmask);

% cells of interest
kk1 = 65;
kk2 = 199; 

mean_nldas_precip = nanmean(nldas.prec,3);

ymax=400;

ms = 30;
figure
plotraster(basin.lonv, basin.latv, mean_nldas_precip, 'Allegheny basin (NLDAS precip mm/day)')
hold on
plot(lon(kk1)-1/16, lat(kk1)-1/16, '.r', 'markersize', ms)
text(lon(kk1), lat(kk1), 'Cell 1', 'fontsize', fs)
plot(lon(kk2)-1/16, lat(kk2)-1/16, '.r', 'markersize', ms)
text(lon(kk2), lat(kk2), 'Cell 2', 'fontsize', fs)
% plot(basin.gage_lon(1), basin.gage_lat(1), '.b', 'markersize', ms)
% text(basin.gage_lon(1), basin.gage_lat(1)+1/16, 'Gage 1', 'fontsize', fs)
% plot(basin.gage_lon(2), basin.gage_lat(2), '.b', 'markersize', ms)
% text(basin.gage_lon(2), basin.gage_lat(2)+1/16, 'Gage 2', 'fontsize', fs)
set(gca, 'fontsize', fs)

figure
subplot(1,2,1)
shadedErrorBar(1:365, squeeze(prior_runoff_ens(kk1,:,:))', {@mean, @std}, ...
    'lineprops', '-r');
grid on
hold on
plot(1:365, true_runoff(kk1,:), 'k-', 'linewidth', lw)
xlim([1,365])
ylim([-1,ymax])
legend('Prior mean (TMPA)','Truth (NLDAS)')
xlabel('Day of year')
ylabel('Runoff (mm/day)')
title('Cell 1')
set(gca, 'fontsize', fs)
subplot(1,2,2)
shadedErrorBar(1:365, squeeze(prior_runoff_ens(kk2,:,:))', {@mean, @std}, ...
    'lineprops', '-r');
grid on
hold on
plot(1:365, true_runoff(kk2,:), 'k-', 'linewidth', lw)
xlim([1,365])
ylim([-1,ymax])
xlabel('Day of year')
ylabel('Runoff (mm/day)')
title('Cell 2')
legend('Prior mean (TMPA)','Truth (NLDAS)')
set(gca, 'fontsize', fs)

%% Alternative: generate Gaussian runoff errors

% This method is RAM-limited bc errors are generated in one batch
M=2000; % number of replicates

% Parameters for generating runoff error covariance
opt.L = 40; % km
opt.T = 5; % days 
opt.RAM = 36;
opt.rho_thres = 0;
opt.progress=1;
opt.alpha = 1; % proportionality constant
!rm ./sc.mat
!rm ./tc.mat
!rm ./rho.mat
opt.SC_savename = './sc.mat';
opt.TC_savename = './tc.mat';
opt.rho_savename = './rho.mat';

% T = 5, 25, 1e6 discharge estimates comparison...

% Generate normally-distributed errors
mu1 = zeros(n*nt,1); % mean of runoff errors

% tmp1 = true_runoff(1:n,1:nt)';
% tmp2 = flipud(tmp1);
% x4cov1 = reshape(tmp2', n*nt, 1);

x4cov = fliparoo(true_runoff(1:n,1:nt)', nt-1); % ok

[sigma1, rho1] = calc_yang_covariance2(basin, opt, n, nt-1, x4cov); % covariance of runoff errors
sqrt(diag(sigma1))-x4cov; % check, should be zeros

figure
imagesc(sigma1(1:3*n, 1:3*n))
% imagesc(sigma1)
title('(Part of) error cov matrix')
colorbar

% tmp = mvnrnd(ones(M,1), eye(M), M);

% tmp = mvnrnd(mu1, rho1, M)'; % what if we don't scale by the runoff values?
tmp = mvnrnd(mu1, sigma1, M)'; % runoff errors, in reverse order

runoff_errors = zeros(nt,n,M);
for i=1:M
    runoff_errors(:,:,i) = flipud(reshape(tmp(:,i), n, nt)');
end
runoff_errors = permute(runoff_errors, [2,1,3]);
    
figure
subplot(1,2,1)
imagesc(true_runoff), colorbar
title('true runoff')
subplot(1,2,2)
imagesc(std(runoff_errors, [], 3)), colorbar
title('Std of runoff errors')

% Generate prior
prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
    prior_runoff_ens(:,:,mm) = true_runoff + runoff_errors(:,:,mm);
%     prior_runoff_ens(:,:,mm) = tmpa_runoff_prior'.*lognrnd(mu1, sigma1, n, nt);
end
disp('Generated prior ensemble')

% Plot prior runoff maps as a check
figure
t=[20,30];
for i=1:2
    prior_runoff_map_ENS = make_map(basin, squeeze(mean(prior_runoff_ens(:,t(i),:),3)));
    subplot(1,2,i)
    plotraster(basin.lonv, basin.latv, prior_runoff_map_ENS, ['Prior ensemble mean (day ' num2str(t(i)) ')'])
%     caxis([0,6])
end
colormap cool

% Check that errors have desired properties using our suite of
% error-checking functions (not fully working/understood)

analyze_errors_func_M(prior_runoff_ens, true_runoff, 0) % add errors

% Check variogram in R (write out data from MATLAB)

distmat = basin.distmat;
mask = basin.mask;
gridlon = basin.grid_lon;
gridlat = basin.grid_lat;
lon = basin.lonv;
lat = basin.latv;
save('./allegheny_data/errors/runoff_errors_for_variogram.mat', 'runoff_errors', 'lon','lat','mask','distmat','gridlon','gridlat')

roe = make_map(basin, runoff_errors(:,1,1));
figure,plotraster(basin.lonv, basin.latv, roe, 'runoff error')

% reshape to 2D matrix for writing to file
[nr,nc,nt] = size(runoff_error);
runoff_error_2d = reshape(runoff_error, nr*nc, nt);
dlmwrite('./ohio_data/nldas_tmpa_error.txt', runoff_error_2d, 'delimiter', '\t')
dlmwrite('./ohio_data/distmat.txt', basin.distmat, 'delimiter', '\t')

for tt=1:365
    [x,y,z] = grid2xyz(basin.lonv', basin.latv', runoff_error(:,:,tt));
    dlmwrite(['./ohio_data/runoff_error_xyz/xyz' num2str(tt) '.txt'], [x,y,z], 'delimiter', '\t')
end



% analyze_errors_func(tmpa_runoff_prior, true_runoff', 0) % add errors
% analyze_errors_func(tmpa_runoff_prior, true_runoff', 1) % mult errors

%% Generate discharge "measurements" (Ensemble)

rng(704753262,'twister')

true_discharge = gage;
swot = 1;
percent_error = 15;

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

lrv = lognrnd(mu1, sigma1, 1e4, 1);
rv = log(lrv);

% figure
% plot(tv, true_discharge)
% hold on
% plot(tv, error_corrupted_discharge_meas)
% legend('True discharge','Synthetic measurements')
% xlabel('Time')
% ylabel('Q (mm/day)')

% Plot measurements for each gage
figure
lw = 2;
ind = 1;
for gg=[2]
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

cal_ind = 1:14;
val_ind = [];

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
basin.mask = flipud(basin.mask); % needs to be upside down after this command

figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be upside down here

opt.gauss = 1; % use additive Gaussian errors (vs. lognormal mult errors)
local = 0;
opt.pbs = 0; % option for particle batch smoother (not working)
opt.multerrors = 0; % using multiplicative (log) update
tic
if local
    % with domain localization
    [post_runoff, w1, w2, maxv, maxv_sub] = ISR_domain_wrapper(prior_runoff_ens, ...
        HH(cal_ind,:,:), error_corrupted_discharge_meas(:,cal_ind), s, basin, Cv, opt, tv);
else
    % without localization
    opt.loc = 0;
    [post_runoff, ~, ~, p] = ISR_EnKF(prior_runoff_ens, ...
        HH(cal_ind,:,:), error_corrupted_discharge_meas(:,cal_ind), s, basin, Cv, opt);
    % One gage only
%     post_runoff = ISR_EnKF(prior_runoff_ens, ...
%             HH(1,:,:), error_corrupted_discharge_meas(:,1), s, C.basin, Cv, opt);    
end
toc

basin.mask=flipud(basin.mask);

truth = struct();
truth.total_runoff = nldas_runoff_true';
truth.true_runoff = nldas_runoff_true';
basin.true_runoff = truth.true_runoff;
mean_post_runoff = mean(post_runoff,3);

%    save('./allegheny_data/final/ens_M200_m1s2L40T5_lognorm_daily_mult_0.mat', 'mean_post_runoff', ...
%        'truth', 'error_corrupted_discharge_meas','basin', 'cal_ind', 'val_ind');

%% How does estimated discharge compare to observations?

prior.mean_runoff = mean(prior_runoff_ens,3);
posterior.mean_runoff = mean(post_runoff,3);

% posterior.mean_runoff(posterior.mean_runoff<0) = 0;

prior.meanQ = state_model_dumb(prior.mean_runoff', HH);
posterior.meanQ = state_model_dumb(posterior.mean_runoff', HH);

mcal=14;
prior.Qens = zeros(nt,mcal,M);
posterior.Qens = zeros(nt,mcal,M);
for kk=1:M
    prior.Qens(:,:,kk) = state_model_dumb(prior_runoff_ens(:,:,kk)', HH(cal_ind,:,:));
    posterior.Qens(:,:,kk) = state_model_dumb(post_runoff(:,:,kk)', HH(cal_ind,:,:));
end

% gi=1:nt;
% gi=30:(nt-(k+1));
gi=(k+2):(nt-(k+1));
gagei = cal_ind(2);

kge_prior = plot_discharge_ts(prior.meanQ(:,gagei), true_discharge(:,gagei))

% Plot means
figure
plot(prior.meanQ(gi,gagei), 'b-*', 'linewidth', lw)
hold on
kge_post = plot_discharge_ts(posterior.meanQ(gi,gagei), true_discharge(gi,gagei))
% title('Discharge (M=200, 10 days, no err)')
plot(error_corrupted_discharge_meas(gi,gagei), 'k.', 'markersize', 20)
% xlim([100,160])
ylim([0,1000])
legend('Prior','Truth','Posterior','Observations')

% Plot ensemble members
figure
h1=plot(squeeze(prior.Qens), 'b-', 'linewidth', 1);
hold on
h2=plot(squeeze(posterior.Qens), 'r-', 'linewidth', 1);
h3=plot(error_corrupted_discharge_meas(:,gagei), 'k.', 'markersize', 20);
h4 = plot(true_discharge(:,gagei), 'k-', 'linewidth', lw);
legend([h1(1), h2(1), h3(1), h4], 'Prior','Posterior','Meas', 'Truth')
title('Prior and posterior discharge ensembles')
xlabel('Days')
ylim([-1000,1000])
ylabel('Discharge (mm/day)')
set(gca, 'fontsize', fs)

% Percent negative
100*sum(posterior.mean_runoff(:)<0)/numel(posterior.mean_runoff) % percent of mean <0
100*sum(post_runoff(:)<0)/numel(post_runoff) % percent of ensemble <0

% Calculate median runoff NSE
priornse = zeros(n,1);
postnse = zeros(n,1);
for kk=1:n
    priornse(kk) = myNSE(true_runoff(kk,gi)', prior.mean_runoff(kk,gi)');
    postnse(kk) = myNSE(true_runoff(kk,gi)', posterior.mean_runoff(kk,gi)');
end
disp(['median NSE: ' num2str(median(priornse))])
disp(['median NSE: ' num2str(median(postnse))])

% Plot mean runoff at a cell
kk=39;
figure
plot(true_runoff(kk,:), 'k','linewidth', lw)
hold on, plot(posterior.mean_runoff(kk,:), 'r-','linewidth', lw)
plot(prior.mean_runoff(kk,:), 'b*-', 'linewidth', lw)
legend('true','post', 'prior')
myNSE(true_runoff(kk,gi)', prior.mean_runoff(kk,gi)')
myNSE(true_runoff(kk,gi)', posterior.mean_runoff(kk,gi)')

% Runoff ensemble at a cell

%% Plot runoff map

nse_min = -0;

figure
subplot(1,2,1)
[priornse, kge, rmse, nsemap1] = plot_gofmaps(basin, mean(prior_runoff_ens,3)', nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', 20)
title('TMPA Prior NSE')
caxis([nse_min,1])
subplot(1,2,2)
[postnse, kge, rmse, nsemap2] = plot_gofmaps(basin, posterior.mean_runoff', nldas_runoff_true', gi);
hold on 
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', 20)
title('Ens PostNSE')
caxis([nse_min,1])

figure
subplot(2,1,1)
histogram(priornse)
xlim([-5,1])
title('Prior NSE')
subplot(2,1,2)
histogram(postnse)
title('Posterior NSE')
xlim([-5,1])

%% Save results

save('./allegheny_data/results/Ens_AG_TMPA_swot15_m0a1L40T5.mat', '-v7.3')
