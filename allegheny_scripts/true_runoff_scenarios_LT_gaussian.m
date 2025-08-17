%% Evaluating the effect of correlation in the true runoff
%
% 8/3/2023 JRS
% Running ISR for multiple simulated runoff scenarios, generated to have
% different correlations
% 
% This time, I am generating Gaussian runoff, instead of lognormal
% That way, it will be consistent w the PW13 and Y20 ISR formulations

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

A=load('./allegheny_data/setup-2-gage.mat');
load('./allegheny_data/setup-233-gage.mat'); % only one gauge is actually used in this analysis
basin.distmat = A.basin.distmat;
clearvars A
load('./allegheny_data/alleg_nldas_nanfilled.mat')
load('./allegheny_data/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

true_discharge = gage;
[nt,m] = size(true_discharge);

%% Plot tmpa prior vs. nldas (true) runoff (Y20 Figure 6)

% TMPA prior, NLDAS truth (run as a check)

figure
t=75:77;
cmax = 4;
for i=1:3
    subplot(2,3,i)
    plotraster(basin.lonv, basin.latv, tmpa.runoff(:,:,t(i)), ['TMPA runoff (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
    subplot(2,3,3+i)
    plotraster(basin.lonv, basin.latv, nldas.runoff(:,:,t(i)), ['NLDAS runoff (day ' num2str(t(i)) ')'])
    caxis([0,cmax])
end
colormap cool

%% 

% Runoff values have been simulated using CoSMoS
% nreps replicates for each LT scenario

% Assemble runoff prior
runoffdir = './allegheny_data/runoff_sim_norm5';
% Lvals = [0.01, 3, 6 ,9, 12, 15]; % grid cells (12.31 km per grid cell)
% Tvals = [0.01, 1, 2, 3, 4]; % days

Lvals = linspace(0.01, 50, 50); % grid cells
Tvals = linspace(0.01, 25, 50); % days

% Lvals range from 0 to Lstar
% Tvals range from 0 to Tstar
nreps = 5;
runoff_sim_all = zeros(n, nt, length(Lvals), length(Tvals), nreps);
discharge_true = zeros(nt,m, length(Lvals), length(Tvals), nreps);
tic
for i=1:length(Lvals)
    for j=1:length(Tvals)
        for repl = 1:nreps
            % read in CoSMoS runoff
          runoffname = fullfile(runoffdir, ['L' num2str(Lvals(i), '%.3f') '_T' num2str(Tvals(j), '%.3f') '_' num2str(repl) '.txt']);
%             runoffname = fullfile(runoffdir, ['L' num2str(round(Lvals(i),13), 15) '_T' num2str(Tvals(j)) '_' num2str(repl,3) '.txt']);
            A = dlmread(runoffname); 
            dim = sqrt(size(A,2));
            B = reshape(A(1:nt,:), nt, dim, dim);
            C = permute(B, [2,3,1]);
            D = C(1:20,:,:); % basin area is rect, not square
            runoff_sim = reshape(D, 20*21, nt);
            basinmask = basin.mask; % no nan values
            basinmask(isnan(basin.mask)) = 0;
            basinmask = logical(basinmask);
            runoff_sim_all(:,:,i,j,repl) = runoff_sim(basinmask, :);
            % calculate true discharge
            discharge_true(:,i,j,repl) = state_model_dumb(runoff_sim(basinmask, :, :)', HH(1,:,:));
        end
    end
end
toc
% took about 8 minutes ^^

discharge_true = zeros(nt,m, length(Lvals), length(Tvals), nreps);
for i=1:length(Lvals)
    for j=1:length(Tvals)
        for repl = 1:nreps
            discharge_true(:,:,i,j,repl) = state_model_dumb(runoff_sim_all(:,:,i,j,repl)', HH);
%             discharge_true(:,i,j,repl) = state_model_dumb(runoff_sim_all(:,:,i,j,repl)', HH(1,:,:));
        end
    end
end

% save('./allegheny_data/simulated_runoff_norm5.mat', 'runoff_sim_all', 'Lvals', 'Tvals', 'nreps','-v7.3')
% save('./allegheny_data/simulated_discharge_norm5.mat', 'discharge_true')
% save('./allegheny_data/simulated_discharge_norm5_m233.mat', 'discharge_true')

% load('./allegheny_data/simulated_runoff_norm5.mat')
% load('./allegheny_data/simulated_discharge_norm5.mat')
%% Assemble runoff prior

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

%% Plot discharge "measurements" (PW13)

lw = 2;
ms = 30;
fs = 16;
figure
plot(tv, discharge_true(:,:,1,50,1), 'linewidth', lw)
hold on 
plot(tv, discharge_true(:,:,1,50,1), '.', 'markersize', ms)
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
cov = 1; % coefficient of variation
R = 0;

runoff_prior = 6*ones(nt,n); % mean true runoff is 6, also

post_runoff = zeros(nt, n, length(Lvals), length(Tvals), nreps);
for i=1:1%length(Lvals)
    for j=1:length(Tvals)
        for repl = 1:nreps
%             post_runoff(:,:,i,j,repl) = ISR_PW13(tmpa_runoff_prior, ...
%                 HH(1,:,:), discharge_true(:,i,j,repl), s, 'proportional', cov, R);
            post_runoff(:,:,i,j,repl) = ISR_PW13(runoff_prior, ...
                HH, discharge_true(:,:,i,j,repl), s, 'const_diag', cov, R);            
        end
    end
end

% save('./allegheny_data/post_runoff_PW13_08072023.mat', 'post_runoff')
% post_runoff = load('./allegheny_data/post_runoff_PW13_08042023.mat')
% post_runoff = post_runoff.post_runoff;

% post_runoff(:,:,i,j,repl) = ISR_PW13(ones(nt,n), ...
%     HH(1,:,:), discharge_true(:,i,j,repl), s, 'proportional', cov, R);
            
%% Show the range of runoff priors

i=50;j=50;

ro = mean(runoff_sim_all(:,:,i,j,:),5);
ro_ts = mean(ro,1);
ro_map = zeros(20,21, nt);
for tt=1:nt
    ro_map(:,:,tt) = make_map(basin, ro(:,tt));
end
ro_mapm = mean(ro_map,3);

% Looking at averages
figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, ro_mapm, 'Mean Runoff')
subplot(1,3,[2,3])
plot(tv, ro_ts, 'linewidth', lw)
title(['L = ' num2str(12.31*Lvals(i),3), ' km, T = ' num2str(Tvals(j)) ' days'])
xlabel('Time')
ylabel('Runoff (mm/day)')
set(gca, 'fontsize', fs)

% Looking at one cell, one time snapshot
figure
subplot(1,3,1)
plotraster(basin.lonv, basin.latv, ro_map(:,:,74), 'Runoff (mm/day) (day 74)')
subplot(1,3,[2,3])
plot(tv, ro(206,:), 'linewidth', lw)
title(['Cell 206, ' 'L = ' num2str(12.31*Lvals(i),3), ' km, T = ' num2str(Tvals(j)) ' days'])
xlabel('Time')
ylabel('Runoff (mm/day)')
set(gca, 'fontsize', fs)

% Based on estimates for Allegheny, Kanawha, and UpTenn,
% realistically, we're likely to see values of L = (50-200 km) 
% and T = (5-15 days), so maybe focus on this region

% What does discharge look like? 

figure
subplot(1,2,1)
plot(tv, mean(runoff_sim_all(:,:,i,j,1),1))
title('Runoff')
subplot(1,2,2)
plot(tv, discharge_true(:,:,i,j,1))
title('Discharge')


%% Look closely at results

i=50;j=1;
ro = mean(runoff_sim_all(:,:,i,j,:),5); % true runoff
post_runoff_repavg = mean(post_runoff(:,:,i,j,:),5)'; % posterior runoff

% Plot runoff ts for a cell
kk=206; % cell number
figure
plot(tv, post_runoff_repavg(kk,:), 'r', 'linewidth', lw)
hold on
plot(tv, ro(kk,:), 'k', 'linewidth', lw)
legend('Posterior', 'Truth')
title(['L = ' num2str(12.31*Lvals(i),3), ' km, T = ' num2str(Tvals(j)) ' days'])
set(gca, 'fontsize', fs)

% NSE histograms
nse = zeros(n,1);
for kk=1:n
    nse(kk) = myNSE(ro(kk,:)', post_runoff_repavg(kk,:)');
end

figure, histogram(nse)
title(['L = ' num2str(12.31*Lvals(i),3), ' km, T = ' num2str(Tvals(j)) ' days'])
xlabel('NSE')
set(gca, 'fontsize', fs)
xlim([-1,1])

% Inversion increment


%% Plot overview of results

truth = struct();
i = 1;
j = 1;

% Average over a particular LT combo
post_runoff_repavg = mean(post_runoff(:,:,i,j,:),5);
truth.total_runoff = mean(runoff_sim_all(:,:,i,j,:),5);

gi = (k+1):nt-(k+1);
basin.true_runoff = truth.total_runoff;
plt_ISR_results_overview(basin, runoff_prior', post_runoff_repavg', truth, tv, gi)

% posterior NSE
nse = zeros(n,1);
post_runoff_nse = zeros(length(Lvals), length(Tvals));
for i=1:length(Lvals)
    for j=1:length(Tvals)
        post_runoff_repavg = mean(post_runoff(:,:,i,j,:),5);
        truth.total_runoff = mean(runoff_sim_all(:,:,i,j,:),5);        
        for kk=1:n
            nse(kk) = myNSE(truth.total_runoff(kk,:)', post_runoff_repavg(:,kk));
        end
        post_runoff_nse(i,j) = mean(nse);
    end
end

figure
imagesc(Tvals, 12.31*Lvals, post_runoff_nse), colorbar
title('Mean NSE')
ylabel('L (km)')
xlabel('T (days)')
set(gca, 'fontsize', fs)
%%
for gridcell=81:150
    post_runoff_repavg = mean(post_runoff(:,:,i,j,:),5);
    truth.total_runoff = mean(runoff_sim_all(:,:,i,j,:),5); 
    figure
    plot(truth.total_runoff(gridcell,:)')
    hold on
    plot(post_runoff_repavg(:,gridcell))
end

% the posterior values are different in each case bc the discharge meas are
% different in each case. The runoff initial guess is irrelevant

%% Figure 1

% i=5; % 50 km
i=5; % 250 km
% j=3; % 1 day
j=3; % 10 days
ro = mean(runoff_sim_all(:,:,i,j,:),5);
ro_ts = mean(ro,1);
ro_map = zeros(20,21, nt);
for tt=1:nt
    ro_map(:,:,tt) = make_map(basin, ro(:,tt));
end
ro_mapm = mean(ro_map,3);

ro_mapm_small = ro_mapm;
ro_ts_small = ro_ts;

% Looking at averages

Tstar = k+1;
Lstar = sqrt(basin.gage_area(1))/1e3;

figure

subplot(2,3,1)
plotraster(basin.lonv, basin.latv, ro_mapm_small, 'Mean Runoff (L = 50 km)')

subplot(2,3,2)
plotraster(basin.lonv, basin.latv, ro_mapm, 'Mean Runoff (L = 250 km)')

subplot(2,3,3)
a = mean(post_runoff_nse,2);
plot([Lstar Lstar],[min(a) 1], 'k', 'linewidth', lw)
hold on
plot(12.31*Lvals, a, 'linewidth', lw)
plot(12.31*Lvals(5), a(5), '.r', 'markersize', ms)
text(12.31*Lvals(5)+10, a(5)-0.025, 'L = 50 km day', 'fontsize', fs)
plot(12.31*Lvals(j), a(j), '.r', 'markersize', ms)
text(12.31*Lvals(j)+10, a(j)-0.025, 'L = 250 km', 'fontsize', fs)
text(Lstar+5, min(a) + 0.025, 'L* = 183 km', 'fontsize', fs)
xlabel('L_{true}')
ylabel('NSE')
grid on
ylim([min(a),max(a)])
set(gca, 'fontsize', fs)

subplot(2,3,4)
plot(tv, ro_ts_small, 'linewidth', lw)
title(['T = 1 day'])
xlabel('Time')
ylabel('Runoff (mm/day)')
set(gca, 'fontsize', fs)

subplot(2,3,5)
plot(tv, ro_ts, 'linewidth', lw)
title(['T = 10 days'])
xlabel('Time')
ylabel('Runoff (mm/day)')
set(gca, 'fontsize', fs)

b = mean(post_runoff_nse,1);
subplot(2,3,6)
plot([Tstar Tstar],[0 1], 'k', 'linewidth', lw)
hold on
plot(Tvals, b, 'linewidth', lw)
plot(Tvals(3), b(3), '.r', 'markersize', ms)
text(Tvals(3)+1, b(3)-0.025, 'T = 1 day', 'fontsize', fs)
plot(Tvals(j), b(j), '.r', 'markersize', ms)
text(Tvals(j)+1, b(j)-0.025, 'T = 10 days', 'fontsize', fs)
text(Tstar+0.5, min(b) + 0.025, 'tc = 4 days', 'fontsize', fs)
xlabel('T_{true}')
ylabel('NSE')
grid on
ylim([min(b),max(b)])
set(gca, 'fontsize', fs)

%% Figure 1, simplified

figure

a = mean(post_runoff_nse,2);
subplot(1,2,1)
plot(12.31*Lvals, a, 'linewidth', lw)
xlabel('L_{true} (km)')
ylabel('NSE')
grid on
title('NSE vs. spatial correlation')
ylim([min(a),max(a)])
xlim([0,12.31*50])
set(gca, 'fontsize', fs)

b = mean(post_runoff_nse,1);
subplot(1,2,2)
plot([Tstar Tstar],[0 1], 'k', 'linewidth', lw)
hold on
plot(Tvals, b, 'linewidth', lw)
text(Tstar+0.5, min(b) + 0.025, 'tc = 4 days', 'fontsize', fs)
xlabel('T_{true} (days)')
ylabel('NSE')
grid on
title('NSE vs. temporal correlation')
ylim([min(b),max(b)])
set(gca, 'fontsize', fs)

%% Can we do better by including correlations in the ISR model a la Y20?

s = 2*(k+1);
optionsfile = './allegheny_data/options_alleg.txt';
rho_thres = exp(-2);

% Doing this for two cases:
% 1. L = 50 km, T = 1 day (i = 5, j = 3)
% 2. L = 250 km, T = 10 days (i = 21, j = 21)

% Trying a range of (L,T) assumptions
%i = 5;
%j = 3;
i = 21;
j = 21;

nL = 50;
nT = 50;
Lvals_assumed = linspace(0.01, 12.31*50, nL); % km
Tvals_assumed = linspace(0.01, 25, nT); % days

runoff_prior = 6*ones(nt,n); 
post_runoff_Y20 = zeros(nt, n, nL, nT, nreps);
for repl = 1:nreps
    for ii=1:nL
        for jj=1:nT
            % reset stored correlation files
            !rm './allegheny_data/SC_99.mat' 
            !rm './allegheny_data/TC_99.mat'
            !rm './allegheny_data/rho_99_99.mat'        
            post_runoff_Y20(:,:,ii,jj, repl) = ISR_Y20(runoff_prior, HH(1,:,:), discharge_true(:,i,j,repl), ...
                s, basin, optionsfile, Lvals_assumed(ii), Tvals_assumed(jj), rho_thres);
        end
        clc 
        disp(ii)
    end
end

save('./allegheny_data/Y20_results_high_true_corr_M5.mat', 'post_runoff_Y20', '-v7.3')

%% Looking at the effect of spatial correlation only

% 9/1/2023

optionsfile = './allegheny_data/options_alleg.txt';
rho_thres = exp(-2);

nL = 50;
nT = 50;
Lvals_assumed = linspace(0.01, 12.31*50, nL); % km
Tvals_assumed = linspace(0.01, 25, nT); % days

% Loop over true L values
post_runoff_Y20_L = zeros(nt, n, nL, 2, nreps);
% Make one data series using L_assumed = 0 and one with L_assumed = L_true
for case_ind=1:2
    for repl=1:nreps
        for ii=1:nL
            if case_ind==2
            % reset stored correlation files
            !rm './allegheny_data/SC_99.mat' 
            !rm './allegheny_data/TC_99.mat'
            !rm './allegheny_data/rho_99_99.mat'        
            post_runoff_Y20_L(:,:,ii,2,repl) = ISR_Y20(runoff_prior, HH(1,:,:), discharge_true(:,ii,1,repl), ...
                s, basin, optionsfile, Lvals_assumed(ii), Tvals_assumed(1), rho_thres);
            elseif case_ind==1
               post_runoff_Y20_L(:,:,ii,1,repl) = ISR_PW13(runoff_prior, ...
                    HH(1,:,:), discharge_true(:,ii,1,repl), s, 'const_diag', cov, R); 
            end
        end
        clc
        disp(repl)
    end
end

nse = zeros(n,2);
post_runoff_nse_L = zeros(50,2);
for case_ind=1:2
    for i=1:50
        post_runoff_repavg = mean(post_runoff_Y20_L(:,:,i,case_ind,:),5);
        truth.total_runoff = mean(runoff_sim_all(:,:,i,1,:),5);        
        for kk=1:n
            nse(kk,case_ind) = myNSE(truth.total_runoff(kk,:)', post_runoff_repavg(:,kk));
        end
        post_runoff_nse_L(i,case_ind) = mean(nse(:,case_ind));
    end
end

figure(16)
clf
subplot(1,2,1)
plot(12.31*Lvals, post_runoff_nse_L, 'linewidth', lw)
xlabel('L_{true} (km)')
ylabel('NSE')
grid on
legend('L=0', 'L = L_{true}')
title('NSE vs. spatial correlation')
ylim([0,1])
xlim([0,12.31*50])
set(gca, 'fontsize', fs)

%% Looking at temporal correlation only 

% 9/1/2023
% In the context of SWOT-like measurements

% Loop over true L values
post_runoff_Y20_T = zeros(nt, n, nT, 2, nreps);
% Make one data series using L_assumed = 0 and one with L_assumed = L_true
for case_ind=1:2
    for repl=1:nreps
        for jj=1:nT
            if case_ind==2
            % reset stored correlation files
            !rm './allegheny_data/SC_99.mat' 
            !rm './allegheny_data/TC_99.mat'
            !rm './allegheny_data/rho_99_99.mat'        
            post_runoff_Y20_T(:,:,jj,2,repl) = ISR_Y20(runoff_prior, HH(1,:,:), discharge_true(:,4,jj,repl), ...
                s, basin, optionsfile, Lvals_assumed(4), Tvals_assumed(jj), rho_thres);
            elseif case_ind==1
               post_runoff_Y20_T(:,:,jj,1,repl) = ISR_PW13(runoff_prior, ...
                    HH(1,:,:), discharge_true(:,4,jj,repl), s, 'const_diag', cov, R); 
            end
        end
        clc
        disp(repl)
    end
end

nse = zeros(n,2);
post_runoff_nse_T = zeros(50,2);
for case_ind=1:2
    for j=1:50
        post_runoff_repavg = mean(post_runoff_Y20_T(:,:,j,case_ind,:),5);
        truth.total_runoff = mean(runoff_sim_all(:,:,1,j,:),5);        
        for kk=1:n
            nse(kk,case_ind) = myNSE(truth.total_runoff(kk,:)', post_runoff_repavg(:,kk));
        end
%         figure
%         plot(mean(post_runoff_repavg,2)), hold on, 
%         plot(mean(truth.total_runoff))
%         legend('est','true')
        post_runoff_nse_T(j,case_ind) = mean(nse(:,case_ind));
    end
end

figure(16)
subplot(1,2,2)
plot(Tvals, post_runoff_nse_T, 'linewidth', lw)
xlabel('T_{true} (km)')
ylabel('NSE')
grid on
legend('T=0','T = T_{true}')
title('NSE vs. temporal correlation')
ylim([-1,0])
xlim([0,25])
set(gca, 'fontsize', fs)

%% Looking at temporal correlation when there are missing measurements

L_ind = 20;

post_runoff_Y20_T_swot = zeros(nt, n, nT, 4, 2);
% post_runoff_Y20_T_swot = zeros(nt, n, nT, 4, nreps);
% Using L_assumed = L_true
freq = [2, 3, 5, 10];
for case_ind=1:4
    for repl=1:2
        for jj=1:nT
            % reset stored correlation files
            !rm './allegheny_data/SC_99.mat' 
            !rm './allegheny_data/TC_99.mat'
            !rm './allegheny_data/rho_99_99.mat'     
            allobs = discharge_true(:,:,L_ind,jj,repl);
            swot_obs = NaN(size(allobs));
            swot_obs(1:freq(case_ind):end,:) = allobs(1:freq(case_ind):end,:);
            post_runoff_Y20_T_swot(:,:,jj,case_ind,repl) = ISR_Y20(runoff_prior, HH, swot_obs, ...
                s, basin, optionsfile, Lvals_assumed(L_ind), Tvals_assumed(jj), rho_thres);
            
%              pr = ISR_PW13(runoff_prior, ...
%                      HH, discharge_true(:,:,1,1,repl), s, 'const_diag', cov, R);             
                 
        end
%         clc
        disp(repl)
    end
end

i=1; % spatial correlation index
j=1; % temporal correlation index
case_ind = 1; % sampling frequency index
allobs = discharge_true(:,:,i,j,repl);
swot_obs = NaN(size(allobs));
swot_obs(1:freq(case_ind):end,:) = allobs(1:freq(case_ind):end,:);
figure
outlet_ind = 29;
plot(allobs(:,29))
hold on
plot(swot_obs(:,29), 'r.', 'markersize', 20)
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat'             
pr = ISR_Y20(runoff_prior, HH, swot_obs, ...
                s, basin, optionsfile, Lvals_assumed(i), Tvals_assumed(j), rho_thres);
pd = state_model_dumb(pr, HH);
gageind=1; % 29=outlet
plot_discharge_ts(pd(:,gageind), discharge_true(:,gageind,L_ind,jj,repl))
hold on
plot(swot_obs(:,gageind), 'r.', 'markersize', 20)
plot_discharge_ts(pd, discharge_true(:,:,L_ind,jj,repl))

figure
plot(runoff_sim_all(1,:,L_ind,jj,repl))
hold on
plot(pr(:,1))
legend('post','true')

% plot_discharge_ts(post_discharge_T_swot(:,29,1,1,1), discharge_true(:,29,1,1,1))

nse = zeros(n,4);
post_runoff_nse_T_swot = zeros(size(post_runoff_Y20_T_swot,3),4);
for case_ind=1:4
    for j=1:size(post_runoff_Y20_T_swot,3)
        post_runoff_repavg = mean(post_runoff_Y20_T_swot(:,:,j,case_ind,:),5);
        truth.total_runoff = mean(runoff_sim_all(:,:,1,j,:),5);        
        for kk=1:n
            nse(kk,case_ind) = myNSE(truth.total_runoff(kk,:)', post_runoff_repavg(:,kk));
        end
%         figure
%         plot(mean(post_runoff_repavg,2)), hold on, 
%         plot(mean(truth.total_runoff))
%         legend('est','true')
        post_runoff_nse_T_swot(j,case_ind) = mean(nse(:,case_ind));
    end
end

% Plot runoff for a grid cell
figure
plot(mean(post_runoff_repavg,2)), hold on, 
plot(mean(truth.total_runoff))
legend('est','true')

post_discharge_T_swot = zeros(nt,m,nT,4,2);
% post_discharge_T_swot = zeros(nt,m,nT,nreps);
for repl=1:2
    for case_ind=1:4
        for j=1:nT
            post_discharge_T_swot(:,:,j,case_ind,repl) = state_model_dumb(post_runoff_Y20_T_swot(:,:,j,case_ind,repl),HH);
        end
    end
end
    
plot_discharge_ts(post_discharge_T_swot(:,29,1,1,1), discharge_true(:,29,1,1,1))

figure(16)
subplot(1,2,2)
plot(Tvals, post_runoff_nse_T_swot, 'linewidth', lw)
xlabel('T_{true} (km)')
ylabel('NSE')
grid on
legend('2', '3 days','5 days','10 days')
title('NSE vs. temporal correlation')
% ylim([0,1])
xlim([0,25]) 
set(gca, 'fontsize', fs)



%%
% post_runoff_infL; % estimated using L = 1e6, T = 1e-4;
% should be quite bad. True runoff i=1, j=1 (no corr)

ltc = load('./allegheny_data/Y20_results_low_true_corr_M5.mat', 'post_runoff_Y20');
ltc = ltc.post_runoff_Y20;

htc = load('./allegheny_data/Y20_results_high_true_corr_M5.mat', 'post_runoff_Y20');
htc = htc.post_runoff_Y20;

% ltc plots look good (need to change i,j to make plots, btw)
post_runoff_Y20 = htc; % assign to htc or ltc for plotting

% save('./allegheny_data/Y20_results_high_true_corr_M1.mat', 'post_runoff_Y20')

%% Plot grid showing NSE vs. L,T

% Average over a particular LT combo
truth.total_runoff = mean(runoff_sim_all(:,:,i,j,:),5);
gi = (k+1):nt-(k+1);
basin.true_runoff = truth.total_runoff;
plt_ISR_results_overview(basin, runoff_prior', post_runoff_Y20(:,:,ii,jj)', truth, tv, gi) % ii, jj subscripts
% plt_ISR_results_overview(basin, runoff_prior', post_runoff_infL', truth, tv, gi)

post_runoff_Y20 = mean(htc.post_runoff_Y20,5);
% post_runoff_Y20 = ltc.post_runoff_Y20(:,:,:,:,5);

% Also plot discharge at outlet
prior_discharge = state_model_dumb(runoff_prior, HH(1,:,:));
post_discharge_Y20 = state_model_dumb(post_runoff_Y20(:,:,ii,jj), HH(1,:,:));

figure
plot(tv, prior_discharge, 'b-', 'linewidth', lw)
hold on
plot(tv, post_discharge_Y20, 'r-', 'linewidth', lw)
plot(tv, discharge_true(:,i,j), 'k.', 'markersize', 15)
xlabel('Time')
ylabel('Discharge (mm/day)')
legend('Prior','Posterior','Truth/Obs')
grid on
set(gca, 'fontsize', fs)

% Find for which L, T assumption did we get the best runoff estimates

nsemat = zeros(nL, nT, n);
ro = mean(runoff_sim_all(:,:,i,j,:),5); % true runoff
for ii=1:nL
    for jj=1:nT
        post_runoff_repavg = mean(post_runoff_Y20(:,:,ii,jj,:),5)'; % posterior runoff
        for kk=1:n
            nsemat(ii,jj,kk) = myNSE(ro(kk,:)', post_runoff_repavg(kk,:)');
        end
    end
    ii
end

figure
plot(tv, mean(ro,1))
hold on
plot(tv, mean(post_runoff_repavg,1))
legend('true','est')

% Is NSE best for the correct guess of L, T?
figure
imagesc(Tvals_assumed, Lvals_assumed, mean(nsemat,3))
colorbar
xlabel('T (days)')
ylabel('L (km)')
set(gca, 'Ydir', 'normal')
set(gca, 'fontsize', fs)
% Not exactly

% Calculate prior NSE (should be zero by definition)
nse = zeros(n,1);
for kk=1:n
    nse(kk) = myNSE(ro(kk,:)', runoff_prior); % prior nse
end

%% Averaging over L, plot NSE vs. T

figure
subplot(1,2,1) % nsemat(L, T)
plot(Tvals, mean(mean(nsemat,3),1), 'r-', 'linewidth', lw) 
% title(['NSE vs. assumed T (T_{true} = ' num2str(Tvals(j)) ' days)'])
hold on
plot([Tvals(j) Tvals(j)],[min(mean(mean(nsemat,3),1)) max(mean(mean(nsemat,3),1))],'k', 'linewidth', lw)
title('NSE vs. assumed T')
xlabel('T (days)')
ylabel('Basin mean NSE')
grid on
set(gca, 'fontsize', fs)

% Averaging over T, plot NSE vs. L
subplot(1,2,2)
plot(12.31*Lvals, mean(mean(nsemat,3),2), 'r-', 'linewidth', lw) 
hold on
plot([12.31*Lvals(i) 12.31*Lvals(i)],[min(mean(mean(nsemat,3),2)) max(mean(mean(nsemat,3),2))],'k', 'linewidth', lw)
title('NSE vs. assumed L')
xlabel('L (km)')
ylabel('Basin mean NSE')
grid on
set(gca, 'fontsize', fs)

%% Look at runoff estimate for a given grid cell

% Get cell number from latlon
cell_lat = 41.5;
cell_lon = -79.25;
[~, lons, lats] = make_map(basin, runoff_prior(1,:));
aa = find(abs(lats - cell_lat)<=res/2);
bb = find(abs(lons - cell_lon)<res);
[cc,dd] = ismember(aa, bb);
kk = aa(cc);
[lons(kk), lats(kk)];
kk = find(aa);

kk = kk(1); % cell number

ii=5;
jj=3;
figure
plot(tv, runoff_prior(:,kk), 'b-', 'linewidth', lw)
hold on
plot(tv, post_runoff_Y20(:,kk,ii,jj), 'r-', 'linewidth', lw)
plot(tv, truth.total_runoff(kk,:), 'k-', 'linewidth', lw)
xlabel('Time')
ylabel('Runoff (mm/day)')
legend('prior','post','true')
nse1 = myNSE(truth.total_runoff(kk,:), post_runoff_Y20(:,kk,ii,jj)');
title(['NSE = ' num2str(nse1)])



