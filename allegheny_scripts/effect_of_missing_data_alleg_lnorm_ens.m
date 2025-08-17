% Determining how missing data affect ISR estimates
%
% 9/3/2023 JRS
% Using lognormal runoff errors and the Ensemble ISR method
%
% Started this, but didn't finish...

%% Load inputs

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

%% Load simulated runoff data

runoffdir = './allegheny_data/runoff_sim_norm5';

Lvals = linspace(0.01, 50, 50); % grid cells
Tvals = linspace(0.01, 25, 50); % days

Tvals = Tvals([1,3,7,11,21,50]);
Lvals = Lvals([1,5,9,17,25,41]);

nT = length(Tvals);
nL = length(Lvals);

nreps=1;
runoff_sim_all = zeros(n, nt, length(Lvals), length(Tvals), nreps);
discharge_true = zeros(nt,m, length(Lvals), length(Tvals), nreps);

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
            discharge_true(:,:,i,j,repl) = state_model_dumb(runoff_sim(basinmask, :, :)', HH);
        end
    end
end
toc

% load('./allegheny_data/simulated_runoff_norm1.mat')
save('./allegheny_data/simulated_runoff_norm1.mat', 'discharge_true','runoff_sim_all', 'Lvals', 'Tvals', 'nreps','-v7.3')

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

%% Y20 ISR parameters

s = 2*(k+1);
cov = 1; % coefficient of variation
R = 0;
optionsfile = './allegheny_data/options_alleg.txt';
rho_thres = exp(-2);
runoff_prior = 6*ones(nt,n); % mean true runoff is 6, also

%% How do the runoff estimates change as fewer data are available?

L_ind = 1; % no corr
T_ind = 1; % no corr

% Choose measurement scenario
freq = [1, 2, 3, 5, 10];
case_ind = 1; 
allobs = discharge_true(:,:,L_ind,T_ind,repl);
swot_obs = NaN(size(allobs));
swot_obs(1:freq(case_ind):end,:) = allobs(1:freq(case_ind):end,:);

% Do Y20 ISR
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat'     
tic
pr = ISR_Y20(runoff_prior, HH, swot_obs, ...
    s, basin, optionsfile, Lvals(L_ind), Tvals(T_ind), rho_thres);
toc
pd = state_model_dumb(pr,HH);

plot_discharge_ts(pd, discharge_true)
for gg=1:m
    kge(gg)=myKGE(discharge_true(32:336,gg), pd(32:336,gg));
end

%% Simulated an AR1 time series

% two cell basin

nreps = 10;
n=2;
xsim = zeros(nt,nreps);
xsim(1,:) = 0;
H = [0 1 1 0 0 0;0,0,0,1,1,0;0,0,0,0,0,1];
rho = 0.5;

basin.mask = [1,1];
basin.mask_gage = [1,1];
basin.flow_length_gage = [10000,0];
flow_vel = 0.1;
timestep = 86400;
[HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, flow_vel, timestep);

% simulate truth
xsim = zeros(nt,n,nreps);
for i=2:nt-1
    xsim(i+1,:,:) = rho*xsim(i,:,:)+randn(1,2,nreps);
end
ysim = zeros(nt,nreps);
for repl=1:nreps
    ysim(:,repl) = state_model_dumb(xsim(:,:,repl), HH);
end
acf(ysim(k+2:end,1),2) % ysim has same corr as xsim 

% meas scenario
freq = [1, 2, 3, 5, 10];
case_ind = 5; 
allobs = ysim(:,1);
swot_obs = NaN(size(allobs));
swot_obs(1:freq(case_ind):end,:) = allobs(1:freq(case_ind):end,:);

% Do Y20 ISR

% Make a plot showing how discharge estimate improves as correlation incr,
% given there is correlation in the runoff data (T)
runoff_prior = ones(nt,n);
yprior = state_model_dumb(runoff_prior, HH);
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat' 
basin.distmat = [0,10000;10000,0];
pr = ISR_Y20(runoff_prior, HH, swot_obs, ...
    s, basin, optionsfile, 10, 10, rho_thres); % L, T
ypred = state_model_dumb(pr,HH);
plot_discharge_ts(ypred, ysim(:,1))
hold on
plot(yprior)
plot(swot_obs, 'r.', 'markersize', 20)
legend('posterior','truth','prior','obs')
plot_discharge_ts(xsim(:,1,1), pr(:,1)) % plot runoff at cell 1
plot_discharge_ts(ysim(:,2,1), pr(:,2)) % plot runoff at cell 2

P = eye(6);
P(13:7:end) = rho;
P(1,5)=rho^2;
P(2,6)=rho^2;
P(5,1)=rho^2;
P(6,2)=rho^2;
P(4,2)=rho;
P(3,1)=rho;
P(5,3)=rho;
P(6,4)=rho;
K = P*H'/(H*P*H')

%% Repeat, but for Allegheny

gagei=29; % outlet
runoff_prior = 6*ones(nt,n);
yprior = state_model_dumb(runoff_prior, HH(gagei,:,:));
repl = 1;

% Choose measurement scenario
freq = [1, 2, 3, 5, 10];
case_ind = 5; 

% True runoff
tr = runoff_sim_all(:,:,L_ind,T_ind)';

% Do Y20 ISR
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat'     
tic
pr = ISR_Y20(runoff_prior, HH(gagei,:,:), swot_obs, ...
    s, basin, optionsfile, 1e7, 1e7, rho_thres); % L, T
ypred = state_model_dumb(pr,HH(gagei,:,:));
kge=plot_discharge_ts(ypred, allobs);
hold on
plot(yprior)
plot(swot_obs, 'r.', 'markersize', 20)
legend('posterior','truth','prior','obs')

kk =11; % cell 1
plot_discharge_ts(tr(:,kk), pr(:,kk)) % plot runoff at cell 1

% Variance of the posterior runoff is very small
% Even though the true runoff has var=1
% The prior runoff has no variance... perhaps that is why
% Could also be a consequence of only having one gauge

%% Figure for paper

% Load true runoff and discharge data 

%% KGE vs. T

L_ind = 4;

gagei=29; % outlet
runoff_prior = 6*ones(nt,n);
yprior = state_model_dumb(runoff_prior, HH(gagei,:,:));
repl = 1;
freq=10;

% Run three times: once for T=0, T=T_true, T=M
M=1e7;
% That=0.01;

kgeT=zeros(nT,3);

for T_ind=1:nT
    That = M;
    allobs = discharge_true(:,gagei,L_ind,T_ind,repl);
    swot_obs = NaN(size(allobs));
    swot_obs(1:freq:end,:) = allobs(1:freq:end,:);    
    !rm './allegheny_data/SC_99.mat' 
    !rm './allegheny_data/TC_99.mat'
    !rm './allegheny_data/rho_99_99.mat'     
    tic
    pr = ISR_Y20(runoff_prior, HH(gagei,:,:), swot_obs, ...
        s, basin, optionsfile, 12.31*Lvals(L_ind), That, rho_thres); % L, T
    ypred = state_model_dumb(pr,HH(gagei,:,:));
    kgeT(T_ind,3)=plot_discharge_ts(ypred, allobs);
    hold on
    plot(yprior)
    plot(swot_obs, 'r.', 'markersize', 20)
    legend('posterior','truth','prior','obs')
end

%% KGE vs. L

T_ind = 4;

% Run three times: once for L=0, L=L_true, L=M
M=1e7;

kgeL=zeros(nT,3);

for L_ind=1:nL
    Lhat = M;
    allobs = discharge_true(:,gagei,L_ind,T_ind,repl);
    swot_obs = NaN(size(allobs));
    swot_obs(1:freq:end,:) = allobs(1:freq:end,:);    
    !rm './allegheny_data/SC_99.mat' 
    !rm './allegheny_data/TC_99.mat'
    !rm './allegheny_data/rho_99_99.mat'     
    tic
    pr = ISR_Y20(runoff_prior, HH(gagei,:,:), swot_obs, ...
        s, basin, optionsfile, Lhat, Tvals(T_ind), rho_thres); % L, T
    ypred = state_model_dumb(pr,HH(gagei,:,:));
    kgeL(L_ind,3)=plot_discharge_ts(ypred, allobs);
    hold on
    plot(yprior)
    plot(swot_obs, 'r.', 'markersize', 20)
    legend('posterior','truth','prior','obs')
end

save('./allegheny_data/effect_of_missing_09032023.mat')

%% Plot

lw=2;
fs=18;

figure
subplot(1,2,1) % nsemat(L, T)
plot(Tvals, kgeT, '-*','linewidth', lw)
legend('T = 0','T=T_{true}','T=M')
title('Discharge KGE at outlet')
xlabel('T (days)')
ylabel('KGE')
grid on
ylim([0,1])
set(gca, 'fontsize', fs)

subplot(1,2,2) % nsemat(L, T)
plot(12.31*Lvals, kgeL, '-*','linewidth', lw)
legend('L = 0','L=L_{true}','L=M')
title('Discharge KGE at outlet')
xlabel('L (km)')
ylabel('KGE')
grid on
ylim([0,1])
set(gca, 'fontsize', fs)

% ^this should probably be done over an ensemble

%% Plot maps of runoff KGE for particular cases

% L,Ltrue = 200km
% Ttrue = 5 days

L_ind = 5;
T_ind = 4;
runoff_true = runoff_sim_all(:,:,L_ind,T_ind)';

% Case 1/2/3: T=0, T=5 days, T=M
case1.allobs = discharge_true(:,gagei,L_ind,T_ind,repl);
case1.swot_obs = NaN(size(case1.allobs));
case1.swot_obs(1:freq:end,:) = case1.allobs(1:freq:end,:);    
That = 0.01;
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat'     
tic
case1.pr = ISR_Y20(runoff_prior, HH(gagei,:,:), case1.swot_obs, ...
    s, basin, optionsfile, 12.31*Lvals(L_ind), That, rho_thres); % L, T
case1.ypred = state_model_dumb(case1.pr,HH(gagei,:,:));

% Case 2
case2.allobs = discharge_true(:,gagei,L_ind,T_ind,repl);
case2.swot_obs = NaN(size(case2.allobs));
case2.swot_obs(1:freq:end,:) = case2.allobs(1:freq:end,:);    
That = Tvals(T_ind);
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat'     
tic
case2.pr = ISR_Y20(runoff_prior, HH(gagei,:,:), case2.swot_obs, ...
    s, basin, optionsfile, 12.31*Lvals(L_ind), That, rho_thres); % L, T
case2.ypred = state_model_dumb(case2.pr,HH(gagei,:,:));

% Case 2
case3.allobs = discharge_true(:,gagei,L_ind,T_ind,repl);
case3.swot_obs = NaN(size(case3.allobs));
case3.swot_obs(1:freq:end,:) = case3.allobs(1:freq:end,:);    
That = M;
!rm './allegheny_data/SC_99.mat' 
!rm './allegheny_data/TC_99.mat'
!rm './allegheny_data/rho_99_99.mat'     
tic
case3.pr = ISR_Y20(runoff_prior, HH(gagei,:,:), case3.swot_obs, ...
    s, basin, optionsfile, 12.31*Lvals(L_ind), That, rho_thres); % L, T
case3.ypred = state_model_dumb(case3.pr,HH(gagei,:,:));

% Calculate GoF
case1.nse_runoff = zeros(n,1);
case2.nse_runoff = zeros(n,1);
case3.nse_runoff = zeros(n,1);
for k=1:n
    case1.nse_runoff(k) = myNSE(runoff_true(:,k), case1.pr(:,k));
    case2.nse_runoff(k) = myNSE(runoff_true(:,k), case2.pr(:,k));
    case3.nse_runoff(k) = myNSE(runoff_true(:,k), case3.pr(:,k));
end
case1.nsemap = make_map(basin, case1.nse_runoff);
case2.nsemap = make_map(basin, case2.nse_runoff);
case3.nsemap = make_map(basin, case3.nse_runoff);

%%
figure

subplot(3,3,1)
plotraster(basin.lonv, basin.latv, case1.nsemap, 'Runoff NSE')
caxis([0,0.5])

subplot(3,3,[2,3])
plot(yprior, 'b', 'linewidth', lw)
hold on
plot_discharge_ts(case1.ypred, case1.allobs);
title('Discharge (T=0)')
plot(case1.swot_obs, 'k.', 'markersize', 20)
legend('Posterior','Truth','Prior','Observations')

subplot(3,3,4)
plotraster(basin.lonv, basin.latv, case2.nsemap, 'Runoff NSE')
caxis([0,0.5])

subplot(3,3,[5,6])
plot(yprior, 'b', 'linewidth', lw)
hold on
plot_discharge_ts(case2.ypred, case2.allobs);
title('Discharge (T=5)')
plot(case2.swot_obs, 'k.', 'markersize', 20)
legend('Posterior','Truth','Prior','Observations')

subplot(3,3,7)
plotraster(basin.lonv, basin.latv, case3.nsemap, 'Runoff NSE')
caxis([0,0.5])

subplot(3,3,[8,9])
plot(yprior, 'b', 'linewidth', lw)
hold on
plot_discharge_ts(case3.ypred, case3.allobs);
title('Discharge (T=M)')
plot(case3.swot_obs, 'k.', 'markersize', 20)
legend('Posterior','Truth','Prior','Observations')

%% Plot true runoff

mean_runoff_map = make_map(basin, mean(runoff_true,1));
figure
plotraster(basin.lonv, basin.latv, mean_runoff_map, 'mean runoff')

% For Ohio, if we can estimate the temporal correlation of the true runoff
% data, then that will likely give us the best estimate of posterior runoff
%
% or, if using big L and T is always best, then we should just choose the
% biggest, but I am pretty sure that big LT is only best for this case
