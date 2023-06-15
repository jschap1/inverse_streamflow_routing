% Hetch Hetchy ISR: parameter sensitivity experiments
%
% Testing sensitivity to L, T, mu, sigma of prior runoff errors
% 5/21/2023

clean
cd /hdd/ISR/02_Basins/HH2/
addpath(genpath('/home/jschap/Documents/Research/Software/VICMATLAB/vicmatlab'))
addpath('/hdd/ISR/ISR/SimpleExample/Functions')
addpath('/hdd/ISR/ISR/buildup/')
addpath('/home/jschap/Documents/MATLAB/NearestSymmetricPositiveDefinite/NearestSymmetricPositiveDefinite')
addpath('/home/jschap/Documents/MATLAB/uimage')

% create inputs using basin_prep_HH2.m
load('./Data/unsorted/isr_setup1.mat')

ms = 20;

basin.mask(isnan(basin.mask)) = 0;
basin.mask = logical(basin.mask);

%% Cut to a particular time period

i1 = find(truth.t == datetime(2006,5,1));
i2 = find(truth.t == datetime(2006,5,31));
timesteps = i1:i2; % originally 5000:7000

% looking at April and May 2007

% timesteps = 5000:7000;
true_discharge = truth.discharge_mm(timesteps,:); % water year 2006
truth.total_runoff = truth.total_runoff(:,timesteps);
tv = truth.t(timesteps);

% cropping to the time period from Oct. 1, 2006 to September 30, 2007

%% EnKF implementation

% Generate ensemble of runoff initial guesses
% (can think of this as ensemble of state at initial time x0)
% We will use the true runoff plus errors
% The errors will be multiplicative lognormal errors with some mu, sigma
% (ultimately, this should be done with CoSMoS with realistic errors)

% s = 2*k+1;
% s = k+1;
s = k+18; % 24 hour window
mean1 = 0.2;
std1 = 1;
mu1 = log((mean1^2)/sqrt(std1+mean1^2));
sigma1 = sqrt(log(std1/(mean1^2)+1));
[n, nt] = size(truth.total_runoff);
m = size(true_discharge,2); % number of gages 
        
% Generate synthetic measurements
% ENKF() assumes homoscedastic, uncorrelated lognormal errors here
% we will use unbiased relative error of 15% (COV=0.15, mu = 1)
mean1 = 1;
% Cv = 0.15;
Cv = (0.15)^2; % variance
mu1 = log((mean1^2)/sqrt(Cv+mean1^2));
sigma1 = sqrt(log(Cv/(mean1^2)+1));
v = lognrnd(mu1, sigma1, 1e5, 1);
% figure
% histogram(v)
mean(v)
std(v)

rng(704753262, 'twister') % set seed to ensure reproducible results
error_corrupted_discharge_meas = true_discharge.*lognrnd(mu1, sigma1, nt, m);

%% Load ensemble of runoff errors

L_vals = [0.01, 1.57, 7, 14, 21, 28, 56, 83, 538];
T_vals = [0.01, 1.5, 1.875, 3, 4.5, 6, 12, 18, 117];
COV = 1;
mu_vals = [0.5, 1, 2];
sigma_vals = mu_vals/COV;

for i=1:length(L_vals)
    for j=1:length(T_vals)
        for kk=1:length(mu_vals)
    
plotflag = 0; % no plotting
            
dirname = ['L', num2str(L_vals(i)), '_T', num2str(T_vals(j)), ...
    '_mu', num2str(mu_vals(kk)), '_sigma', num2str(sigma_vals(kk))];
wkdir_orig = '/hdd/ISR/02_Basins/HH2/Data/runoff_priors/';
wkdir = fullfile(wkdir_orig, dirname);

A = load(fullfile(wkdir, 'simulated_runoff_errors_500.mat'));

M = 500; % ensemble size

A.alldata = A.alldata(:,:,1:nt,:);

runoff_errors = reshape(A.alldata, 6*7, nt, M);
runoff_errors = runoff_errors(basin.mask, :, :);

prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
    prior_runoff_ens(:,:,mm) = truth.total_runoff.*runoff_errors(:,:,mm);
end

% Plot prior ensemble
lw = 2;
tmp = squeeze(mean(prior_runoff_ens,1));

if plotflag
    figure
    h1 = plot(tv, tmp, 'cyan');
    hold on
    h2 = plot(tv, mean(tmp,2), '--*b', 'linewidth', lw);
    h3 = plot(tv, mean(truth.total_runoff,1), 'k', 'linewidth', lw);
    xlabel('Time')
    ylabel('Runoff (mm/hour)')
    legend([h1(1);h2;h3], {'Ensemble','Ens mean','Truth'})
    title('Basin average runoff')
    xlim([datetime(2006,5,13), datetime(2006,5,16)])
    set(gca, 'fontsize', fs)
end

% Run EnKF-ISR (RUN)

% Postulate the discharge measurement error covariance matrix

s=k+1;
tic
options.loc = 0;
options.locL = 5; % cells
options.locT = 72; % hr
options.rho_yy=0;
basin.distmat = calc_distance_matrix(basin.grid_lon(basin.mask), basin.grid_lat(basin.mask));
basin.true_runoff = truth.total_runoff;
basin.t = tv;
basin.true_discharge = true_discharge;

[post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, error_corrupted_discharge_meas, s, basin, Cv, options);
% [post_runoff, w1, w2, maxv, maxv_sub] = ISR_EnKF_domain_wrapper_d2(prior_runoff_ens, ...
%     HH, error_corrupted_discharge_meas, s, basin, Cv, options, tv);

mean_prior_runoff = mean(prior_runoff_ens,3);       
sd_prior_runoff = std(prior_runoff_ens, [], 3);
mean_posterior_runoff = mean(post_runoff, 3);
sd_posterior_runoff = std(post_runoff, [], 3);

saveflag = 0;

fs = 18;
lw = 2;

% confine analysis to a single assimilation window
i1=w1(1);
i2=w2(1);
% xlim([tv(i1), tv(i2)])
% don't look at the spin-up/cooldown part of the window

% average bias of the prior
mean(prior_runoff_ens(:))/mean(truth.total_runoff(:));

% average bias of the posterior
mean(post_runoff(:))/mean(truth.total_runoff(:));

% seems to be a -4.24% bias in the posterior runoff. Why? Shouldn't be
% any... the method seems to introduce it. Can we get rid of it?

% Make maps
% outputs posterior GoF 
%   NEED TO ACCOUNT FOR SPIN-UP/INCOMPLETE UPDATE (use gi arg)

gi = 1+s:nt-s; % spin-up/cool-down handling
[gof] = plot_ensemble_runoff_maps(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, gi, plotflag);
% plot_ensemble_runoff_maps(basin, B.mean_prior_runoff, B.mean_posterior_runoff, truth, tv)

% Save
save(fullfile(wkdir, 'ISR_outputs_1gage.mat'),...
    'tv', 'basin','k','s','truth','gof', 'mean_prior_runoff',...
    'mean_posterior_runoff', 'sd_prior_runoff', 'sd_posterior_runoff');

        end
    end
end

%% Read in results for each case

outputs = table('Size',[243,5],'variableTypes',{'double','double','double','double','double'},...
    'variableNames',{'L','T','mu','wildcard','wildcard2'});

ind = 1;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        for kk=1:length(mu_vals)
    
            dirname = ['L', num2str(L_vals(i)), '_T', num2str(T_vals(j)), ...
                '_mu', num2str(mu_vals(kk)), '_sigma', num2str(sigma_vals(kk))];
            wkdir_orig = '/hdd/ISR/02_Basins/HH2/Data/runoff_priors/';
            wkdir = fullfile(wkdir_orig, dirname);

            outputs.L(ind) = L_vals(i);
            outputs.T(ind) = T_vals(j);
            outputs.mu(ind) = mu_vals(kk);
            
            A = load(fullfile(wkdir, 'ISR_outputs_1gage.mat'));
%             outputs.meanNSE(ind) = mean(A.gof.post.nse);
                        
            outputs.wildcard(ind) = mean(A.gof.post.nse) - mean(A.gof.prior.nse);
            outputs.wildcard2(ind) = mean(A.gof.post.nse);
            
            % Reduction in bias
            prior_bias = mean(mean(A.mean_prior_runoff(:,gi))) - mean(mean(A.truth.total_runoff(:,gi)));
            post_bias = mean(mean(A.mean_posterior_runoff(:,gi))) - mean(mean(A.truth.total_runoff(:,gi)));
            outputs.wildcard(ind) = abs(post_bias - prior_bias);
            
            % Reduction in stddev
%             prior_spread = mean(mean(A.sd_prior_runoff(:,gi)));
%             post_spread = mean(mean(A.sd_posterior_runoff(:,gi)));
%             outputs.wildcard2(ind) = prior_spread - post_spread;            

          % Reduction in COV
            prior_COV = A.sd_prior_runoff(:,gi)./A.mean_prior_runoff(:,gi);
            post_COV = A.sd_posterior_runoff(:,gi)./A.mean_posterior_runoff(:,gi);
            mpriorcv(ind) = mean(prior_COV(:));
            mpostcv(ind) = mean(post_COV(:));
            outputs.wildcard2(ind) = mpriorcv(ind) - mpostcv(ind);            
            
            ind = ind + 1;
            
        end
    end
end

%% Make sensitivity figure for paper

% metrics of interest
% posterior NSE
% change in NSE from prior to posterior
% reduction in uncertainty from prior to posterior
% Prior vs. posterior bias
% Prior vs. posterior standard deviation

% bias ratio = estimate/truth
% relative variability = COV_estimate/COV_true

ms = 20;
fs = 14;

[LL,TT] = meshgrid(L_vals,T_vals);

tc = k+1;
Lstar = sqrt(basin.gage_area/1e6);

LL = LL/Lstar;
TT = TT/tc;

figure

kk=1;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        ind = find(outputs.mu==mu_vals(kk) & outputs.L == L_vals(i) & outputs.T == T_vals(j));
%         Z(i,j) = outputs.medianNSE(ind);
        Z(i,j) = outputs.wildcard(ind);
    end
end
h1 = plot3(LL, TT, Z, '.b', 'markersize', ms);
hold on

kk=2;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        ind = find(outputs.mu==mu_vals(kk) & outputs.L == L_vals(i) & outputs.T == T_vals(j));
%         Z(i,j) = outputs.medianNSE(ind);
        Z(i,j) = outputs.wildcard(ind);
    end
end
h2 = plot3(LL, TT, Z, '.r', 'markersize', ms);

kk=3;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        ind = find(outputs.mu==mu_vals(kk) & outputs.L == L_vals(i) & outputs.T == T_vals(j));
%         Z(i,j) = outputs.medianNSE(ind);
        Z(i,j) = outputs.wildcard(ind);
    end
end
% h3 = surf(LL, TT, Z, '.k'); 
% imagesc(Z)
h3 = plot3(LL, TT, Z, '.k', 'markersize', ms);

grid on
title('Sensitivity to prior parameters')
xlabel('L (km)')
ylabel('T (hours)')

xlabel('L/L^*')
ylabel('T/T^*')

zlabel('Reduction in Bias')
% zlabel('\DeltaNSE')
% zlim([0,1])
legend([h1(1),h2(1),h3(1)], '\mu = 0.5','\mu = 1', '\mu = 2')
set(gca, 'fontsize', fs)

xlim([0,3])
ylim([0,3])

%% Sensitivity figure version 2 (image version)

figure

subplot(1,3,1)
kk=1;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        ind = find(outputs.mu==mu_vals(kk) & outputs.L == L_vals(i) & outputs.T == T_vals(j));
        Z(i,j) = outputs.wildcard(ind);
    end
end
uimagesc(L_vals(1:end-1)/Lstar, T_vals(1:end-1)/tc, Z(1:end-1,1:end-1))
colorbar
xlabel('L/L^*')
ylabel('T/T^*')
title('NSE change (\mu=0.5)')
set(gca, 'fontsize', fs)
% caxis([0,1])

subplot(1,3,2)
kk=2;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        ind = find(outputs.mu==mu_vals(kk) & outputs.L == L_vals(i) & outputs.T == T_vals(j));
        Z(i,j) = outputs.wildcard(ind);
    end
end
uimagesc(L_vals(1:end-1)/Lstar, T_vals(1:end-1)/tc, Z(1:end-1,1:end-1))
colorbar
xlabel('L/L^*')
ylabel('T/T^*')
title('NSE change (\mu=1)')
set(gca, 'fontsize', fs)
% caxis([0,1])

subplot(1,3,3)
kk=3;
for i=1:length(L_vals)
    for j=1:length(T_vals)
        ind = find(outputs.mu==mu_vals(kk) & outputs.L == L_vals(i) & outputs.T == T_vals(j));
        Z(i,j) = outputs.wildcard(ind);
    end
end
uimagesc(L_vals(1:end-1)/Lstar, T_vals(1:end-1)/tc, Z(1:end-1,1:end-1))
colorbar
xlabel('L/L^*')
ylabel('T/T^*')
title('NSE change (\mu=2)')
set(gca, 'fontsize', fs)
% caxis([0,1])

%% Questions to answer

% Does assmilation window make a differnetce?
% In particular, if T is large, and window is small, do we lose some of the
% update? Or does the overlapping windows solve this?

% 1. Figure out if assimilation window makes a difference
% 2. Explain what we see in the 3D plots
% 3. Re-run ISR but store variance metrics (no images)
% 4. Re-generate 3D plots but for change in NSE/KGE, bias ratio
% 5. Make plots for change in variance/COV/relative variability

% 6. Write up for the paper with practical interpretation
% 7. Plan exactly how we will do the full-scale ISR Ohio runs
% 8. Generate sens. figure for case with more gauges (with domain ISR)
% 9. Prove domain ISR is (or isn't) equivalent to regular case under all scenarios 
% 10. Also look at sensitivity to COV. I think with bigger COV we'll get
% better results, since the ensemble will encompass the truth

%% Sensitivity to assimilation window

% Small 
% wkdir = '/hdd/ISR/02_Basins/HH2/Data/runoff_priors/L14_T3_mu2_sigma2';

% Large T
wkdir = '/hdd/ISR/02_Basins/HH2/Data/runoff_priors/L14_T117_mu2_sigma2';

% Small window
% s = k+1;

% Large window
s = 10*(k+1);

gi = 1+s:nt-s; % spin-up/cool-down handling

A = load(fullfile(wkdir, 'simulated_runoff_errors_500.mat'));

A.alldata = A.alldata(:,:,1:nt,:);

runoff_errors = reshape(A.alldata, 6*7, nt, M);
runoff_errors = runoff_errors(basin.mask, :, :);

prior_runoff_ens = zeros(n, nt, M);
for mm=1:M
    prior_runoff_ens(:,:,mm) = truth.total_runoff.*runoff_errors(:,:,mm);
end

% Plot prior ensemble
lw = 2;
tmp = squeeze(mean(prior_runoff_ens,1));

figure
h1 = plot(tv, tmp, 'cyan');
hold on
h2 = plot(tv, mean(tmp,2), '--*b', 'linewidth', lw);
h3 = plot(tv, mean(truth.total_runoff,1), 'k', 'linewidth', lw);
xlabel('Time')
ylabel('Runoff (mm/hour)')
legend([h1(1);h2;h3], {'Ensemble','Ens mean','Truth'})
title('Basin average runoff')
xlim([datetime(2006,5,13), datetime(2006,5,16)])
set(gca, 'fontsize', fs)

% Run EnKF-ISR (RUN)

% Postulate the discharge measurement error covariance matrix

tic
basin.distmat = calc_distance_matrix(basin.grid_lon(basin.mask), basin.grid_lat(basin.mask));
basin.true_runoff = truth.total_runoff;
basin.t = tv;
basin.true_discharge = true_discharge;

[post_runoff, w1, w2] = ISR_EnKF(prior_runoff_ens, HH, error_corrupted_discharge_meas, s, basin, Cv, options);
% [post_runoff, w1, w2, maxv, maxv_sub] = ISR_EnKF_domain_wrapper_d2(prior_runoff_ens, ...
%     HH, error_corrupted_discharge_meas, s, basin, Cv, options, tv);

mean_prior_runoff = mean(prior_runoff_ens,3);       
sd_prior_runoff = std(prior_runoff_ens, [], 3);
mean_posterior_runoff = mean(post_runoff, 3);
sd_posterior_runoff = std(post_runoff, [], 3);

% average bias of the prior
mean(prior_runoff_ens(:))/mean(truth.total_runoff(:))

% average bias of the posterior
mean(post_runoff(:))/mean(truth.total_runoff(:))

% seems to be a -4.24% bias in the posterior runoff. Why? Shouldn't be
% any... the method seems to introduce it. Can we get rid of it?

% Make maps
% outputs posterior GoF
[gof] = plot_ensemble_runoff_maps(basin, mean_prior_runoff, mean_posterior_runoff, truth, tv, gi);
