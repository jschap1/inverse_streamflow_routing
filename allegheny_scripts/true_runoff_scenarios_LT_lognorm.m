%% Evaluating the effect of correlation in the true runoff
%
% 8/1/2023 JRS
% Running ISR for multiple simulated runoff scenarios, generated to have
% different correlations

clear, clc, close all
cd /Volumes/HD3/ISR/inverse_streamflow_routing
addpath(genpath('./src/'))

load('./allegheny_data/setup-2-gage.mat'); % only one gauge is actually used in this analysis
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
% 100 replicates for each LT scenario

% Assemble runoff prior
runoffdir = './allegheny_data/runoff_sim';
% Lvals = [0.01, 3, 6 ,9, 12, 15]; % grid cells (12.31 km per grid cell)
% Tvals = [0.01, 1, 2, 3, 4]; % days

Lvals = linspace(0.01, 50, 50); % grid cells
Tvals = linspace(0.01, 25, 50); % days

% Lvals range from 0 to Lstar
% Tvals range from 0 to Tstar
nreps = 1;
runoff_sim_all = zeros(n, nt, length(Lvals), length(Tvals), nreps);
discharge_true = zeros(nt, length(Lvals), length(Tvals), nreps);
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
            runoff_sim_all(:,:,i,j,repl) = runoff_sim(basinmask, :, :);
            % calculate true discharge
            discharge_true(:,i,j,repl) = state_model_dumb(runoff_sim(basinmask, :, :)', HH(1,:,:));
        end
    end
end
toc
% took about 8 minutes ^^

save('./allegheny_data/runoff_sim/simulated_runoff_d3.mat', 'runoff_sim_all', 'discharge_true', 'Lvals', 'Tvals', 'nreps')

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
plot(tv, discharge_true(:,6,10,1), 'linewidth', lw)
hold on 
plot(tv, discharge_true(:,6,10,1), '.', 'markersize', ms)
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

runoff_prior = ones(nt,n); % mean true runoff is 1, also

post_runoff = zeros(nt, n, length(Lvals), length(Tvals), nreps);
for i=1:length(Lvals)
    for j=1:length(Tvals)
        for repl = 1:nreps
%             post_runoff(:,:,i,j,repl) = ISR_PW13(tmpa_runoff_prior, ...
%                 HH(1,:,:), discharge_true(:,i,j,repl), s, 'proportional', cov, R);
            post_runoff(:,:,i,j,repl) = ISR_PW13(runoff_prior, ...
                HH(1,:,:), discharge_true(:,i,j,repl), s, 'const_diag', cov, R);            
        end
    end
end

save('./allegheny_data/post_runoff_PW13_08072023.mat', 'post_runoff')
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

%% Look closely at results

i=5;j=5;
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
i = 50;
j = 50;

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

% the posterior values are different in each case bc the discharge meas are
% different in each case. The runoff initial guess is irrelevant

%% Figure 1

% i=5; % 50 km
i=21; % 250 km
% j=3; % 1 day
j=21; % 10 days
ro = mean(runoff_sim_all(:,:,i,j,:),5);
ro_ts = mean(ro,1);
ro_map = zeros(20,21, nt);
for tt=1:nt
    ro_map(:,:,tt) = make_map(basin, ro(:,tt));
end
ro_mapm = mean(ro_map,3);

% ro_mapm_small = ro_mapm;
% ro_ts_small = ro_ts;

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

%% Can we do better by including correlations in the ISR model a la Y20?

s = k+1;
optionsfile = './allegheny_data/options_alleg.txt';
rho_thres = exp(-2);

% Doing this for two cases:
% 1. L = 50 km, T = 1 day (i = 5, j = 3)
% 2. L = 250 km, T = 10 days (i = 21, j = 21)

% Trying a range of (L,T) assumptions
% i = 5;
% j = 3;
% i = 21;
% j = 21;

% i = 1;
% j = 1;
% i = 50;
% j = 50;

nL = 5;
nT = 5;
Lvals_assumed = linspace(0.01, 3.57e6, nL); % km
Tvals_assumed = linspace(0.01, 78, nT); % days

runoff_prior = ones(nt,n); 
post_runoff_Y20 = zeros(nt, n, nL, nT);
for ii=1:nL
    for jj=1:nT
        % reset stored correlation files
        !rm './allegheny_data/SC_133.mat' 
        !rm './allegheny_data/TC_15.8.mat'
        !rm './allegheny_data/rho_133_15.8.mat'        
        post_runoff_Y20(:,:,ii,jj) = ISR_Y20(runoff_prior, HH(1,:,:), discharge_true(:,i,j,repl), ...
            s, basin, optionsfile, Lvals_assumed(ii), Tvals_assumed(jj), rho_thres);
    end
    disp(ii)
end

post_runoff_infL; % estimated using L = 1e6, T = 1e-4;
% should be quite bad. True runoff i=1, j=1 (no corr)

ltc = load('./allegheny_data/Y20_results_low_true_corr_M1.mat', 'post_runoff_Y20');
ltc = ltc.post_runoff_Y20;

htc = load('./allegheny_data/Y20_results_high_true_corr_M1.mat', 'post_runoff_Y20');
htc = htc.post_runoff_Y20;

% ltc plots look good (need to change i,j to make plots, btw)
post_runoff_Y20 = htc; % assign to htc or ltc for plotting

% save('./allegheny_data/Y20_results_high_true_corr_M1.mat', 'post_runoff_Y20')

% Average over a particular LT combo
truth.total_runoff = mean(runoff_sim_all(:,:,i,j,:),5);
gi = (k+1):nt-(k+1);
basin.true_runoff = truth.total_runoff;
plt_ISR_results_overview(basin, runoff_prior', post_runoff_Y20(:,:,ii,jj)', truth, tv, gi) % ii, jj subscripts
plt_ISR_results_overview(basin, runoff_prior', post_runoff_infL', truth, tv, gi)

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

% Averaging over L, plot NSE vs. T
figure
subplot(1,2,1) % nsemat(L, T)
plot(Tvals, mean(mean(nsemat,3),1), 'r-', 'linewidth', lw) 
title('NSE vs. assumed T (T_{true} = 10 days)')
xlabel('T (days)')
ylabel('Basin mean NSE')
grid on
set(gca, 'fontsize', fs)

% Averaging over T, plot NSE vs. L
subplot(1,2,2)
plot(12.31*Lvals, mean(mean(nsemat,3),2), 'r-', 'linewidth', lw) 
title('NSE vs. assumed L (L_{true} = 250km)')
xlabel('L (km)')
ylabel('Basin mean NSE')
grid on
set(gca, 'fontsize', fs)



