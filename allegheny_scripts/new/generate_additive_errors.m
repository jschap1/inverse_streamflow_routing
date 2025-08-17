%% Generate additive errors for Ensemble ISR
%
% 10/3/2023 JRS
% This script is meant to be submitted to a remote computer to produce outputs

clear, clc, close all
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
cd /home/jschap/Dropbox/ISR/inverse_streamflow_routing/
addpath(genpath('./src/'))

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
tv = tv(i1:i2);
gage = gage(i1:i2,:);
nldas_runoff = nldas_runoff(:,:,i1:i2);
tmpa_runoff = tmpa_runoff(:,:,i1:i2);
true_runoff = true_runoff(:,i1:i2); 
[nt,m] = size(gage);
tmpa.runoff = tmpa.runoff(:, :, i1:i2);
nldas.runoff  = nldas.runoff(:, :, i1:i2);

%% Load true runoff

% basin.mask = flipud(basin.mask);
figure,plotraster(basin.lonv, basin.latv, basin.mask, 'mask') % should be right side up
saveas(gca, './allegheny_data/figures/basinmask.png')

basin_mask_linear = basin.mask(:);
basin_mask_linear(isnan(basin_mask_linear)) = 0;

tmpa_runoff_linear = reshape(tmpa.runoff, length(basin_mask_linear), nt);
tmpa_runoff_prior = tmpa_runoff_linear(logical(basin_mask_linear),:)';
tmpa_runoff_prior(isnan(tmpa_runoff_prior)) = 0;

nldas_runoff_linear = reshape(nldas.runoff, length(basin_mask_linear), nt);
nldas_runoff_true = nldas_runoff_linear(logical(basin_mask_linear),:)';
nldas_runoff_true(isnan(nldas_runoff_true)) = 0;

%% Generate Gaussian runoff errors

% This method is RAM-limited bc errors are generated in one batch
M=500; % number of replicates

% Parameters for generating runoff error covariance
opt.L = 40; % km
opt.T = 5; % days 
opt.RAM = 14;
opt.rho_thres = 0;
opt.progress=1;
opt.alpha = 1; % proportionality constant
!rm ./sc.mat
!rm ./tc.mat
!rm ./rho.mat
opt.SC_savename = './sc.mat';
opt.TC_savename = './tc.mat';
opt.rho_savename = './rho.mat';

% Generate normally-distributed errors
mu1 = zeros(n*nt,1); % mean of runoff errors

% uniform prior case
totmean = mean(nldas.runoff_mean);
runoff_init = ones(nt, n)*totmean;
x4cov = fliparoo(runoff_init, nt-1); 

% x4cov = fliparoo(tmpa_runoff_prior, nt-1); % proportional to prior mean
% x4cov = fliparoo(true_runoff(1:n,1:nt)', nt-1); % proportional to truth

% covariance of runoff errors
[sigma1, rho1] = calc_yang_covariance2(basin, opt, n, nt-1, x4cov);

% check, should be zeros
diagdiff = sqrt(diag(sigma1))-x4cov;
if abs(max(diagdiff))>0.1
    error('Problem with simulated runoff errors')
end

% or we can assume constant variance, rather than proportional to
% true/prior runoff
% val = 4*mean(x4cov);
% sigma1 = val*rho1;

figure
imagesc(sigma1(1:3*n, 1:3*n))
title('(Part of) error cov matrix')
colorbar
set(gca, 'fontsize', fs)
saveas(gca, './allegheny_data/figures/error_covmat.png')

rng(704753262, 'twister')
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
% saveas(gca, './allegheny_data/figures/runoff_err_std_check.png')

% Generate prior
prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
%     prior_runoff_ens(:,:,mm) = true_runoff + runoff_errors(:,:,mm);
    prior_runoff_ens(:,:,mm) = tmpa_runoff_prior' + runoff_errors(:,:,mm);
%     prior_runoff_ens(:,:,mm) = tmpa_runoff_prior'.*lognrnd(mu1, sigma1, n, nt);
end
disp('Generated prior ensemble (additive errors)')

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
% saveas(gca, './allegheny_data/figures/mean_prior_runoff_snapshots.png')

% Check that errors have desired properties using our suite of
% error-checking functions
% close all
% opt.addmult=0;
% opt.T = 5;
% opt.corr = rho1;
% analyze_errors_func_M(prior_runoff_ens, true_runoff, opt) % add errors
% close all

%% Save prior runoff for Y20

% Y20_prior = prior_runoff_ens(:,:,1); % just one realization
% save('./allegheny_data/errors/m0a1L40T5/Y20_AG_prior_TMPA_first_realization.mat', 'Y20_prior')

%% Save prior runoff ensemble for EnsISR

save('./allegheny_data/errors/m0a1L40T5_prop_prior/ens_AG_m0a1L40T5_uniform_ensemble.mat', ...
    'runoff_errors', 'prior_runoff_ens');

% END runoff error generation

%% Generate (and save) EnsISR prior

% one realization, perturbed with our error model
mu1 = Y20_prior(:);
tmp = mvnrnd(mu1, sigma1, M)'; % runoff errors, in reverse order
runoff_errors_ens = zeros(nt,n,M);
for i=1:M
    runoff_errors_ens(:,:,i) = flipud(reshape(tmp(:,i), n, nt)');
end
runoff_errors_ens = permute(runoff_errors_ens, [2,1,3]);
% Generate prior
Ens_prior = zeros(n, nt, M);
for mm=1:M
    Ens_prior(:,:,mm) = Y20_prior + runoff_errors(:,:,mm);
    prior_runoff_ens(:,:,mm) = tmpa_runoff_prior'.*lognrnd(mu1, sigma1, n, nt);
end
disp('Generated prior ensemble')
save('./allegheny_data/errors/m0a1L40T5/Ens_prior_AGflat.mat', 'Ens_prior')
save('./allegheny_data/errors/m0a1L40T5/Ens_errors_AGflat.mat', 'runoff_errors')

%% Check variogram in R (write out data from MATLAB)

distmat = basin.distmat;
mask = basin.mask;
gridlon = basin.grid_lon;
gridlat = basin.grid_lat;
lon = basin.lonv;
lat = basin.latv;
save('./allegheny_data/errors/m0a1L40T5/runoff_errors_for_variogram_AGflat.mat', 'runoff_errors', 'lon','lat','mask','distmat','gridlon','gridlat')
