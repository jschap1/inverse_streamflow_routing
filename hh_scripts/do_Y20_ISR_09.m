%% Y20 ISR for Hetch Hetchy basin
%
% 6/15/2023 JRS
% Debugging and error checking the method for transfer to Ohio basin
%
% Sensitivity test of Y20 method w.r.t. (L,T)
% Repeat this analysis with ensemble ISR

clear, clc, close all
cd('/hdd/ISR/inverse_streamflow_routing')
addpath(genpath('./src/'))

load('./hh_data/isr_setup_1.mat')
basin.mask = flipud(basin.mask);
% 
% basin.lon = basin.grid_lon(basin.mask);
% basin.lat = basin.grid_lat(basin.mask);
% basin.distmat = calc_distance_matrix(basin.lon, basin.lat);
% basin.true_runoff = truth.total_runoff';

ms = 20;
figure
plotraster(basin.lonv, basin.latv, flipud(basin.flow_accum), 'HH2')
hold on
plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', ms)

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

[nt,m] = size(true_discharge);

%% Generate discharge "measurements"

s = 2*k+1;
mn = 1;
sd = 0.15;
error_corrupted_discharge_meas = generate_discharge_measurements(true_discharge, 'norm', mn, sd, 1,0);
error_corrupted_discharge_meas_swot = generate_discharge_measurements(true_discharge, 'norm', mn, sd, 1,1);

%% Y20 method

% characteristic length and time scales of basin
Lstar = sqrt(basin.gage_area(1)/1000^2);
k = 5;
Tstar = k+1;
rho_thres = exp(-2);
% runoff_init = ones(n,nt);
% runoff_init = 2*truth.total_runoff;
err = zeros(n,nt);
for kk=1:n
    err(kk,:) = mvnrnd(zeros(1,nt), 3*eye(nt));
end
runoff_init = truth.total_runoff + err;
optionsfile = './hh_data/options_hh2.txt';

%%
nL = 11;
nT = 10;
Lfrac = linspace(0.01, 3, nL);
Tfrac = linspace(0.01, 3, nT);
nsemin = zeros(nL,nT);
nsemean = zeros(nL,nT);

tic
for ii = 1:nL
    for jj = 1:nT

        L = Lfrac(ii)*Lstar; % km
        T = Tfrac(jj)*Tstar; % hours
        [runoff, K] = ISR_Y20(runoff_init', HH, error_corrupted_discharge_meas, ...
            s, basin, optionsfile, L, T, rho_thres); % Y20 method

        !rm ./hh_data/SC*
        !rm ./hh_data/TC*
        !rm ./hh_data/rho*

        n = size(runoff,2);
        nse = zeros(n,1);
        for kk=1:n
%             nse(kk) = myKGE(truth.total_runoff(kk,:)', runoff(:,kk));
            nse(kk) = myNSE(truth.total_runoff(kk,:)', runoff(:,kk));
        end

        nsemean(ii,jj) = mean(nse);
        nsemin(ii,jj) = min(nse);

    end
end
toc

%% It looks like larger T is bad, past a certain point

[LL, TT] = meshgrid(Lfrac, Tfrac);

figure
subplot(1,2,1)
surf(LL, TT, nsemean')
colorbar
title('mean nse')
xlabel('L/L*')
ylabel('T/T*')
zlabel('mean nse')
% caxis([-1,1])
% set(gca,'YDir','reverse') 

subplot(1,2,2)
surf(LL, TT, nsemin')
colorbar
% caxis([-1,1])
title('min nse')
xlabel('L/L*')
ylabel('T/T*')
zlabel('min nse')

%% Look at outputs

basin.true_runoff = truth.total_runoff;
gi = (k+1):nt-(k+1);
% basin.mask = flipud(basin.mask);
plt_ISR_results_overview(basin, runoff_init, runoff', truth, tv, gi)

%% Plot runoff at each grid cell

figure
lw = 2;
for k=1:n
    subplot(5,4,k)
    plot(tv, runoff_init(k,:), 'b', 'linewidth', lw)
    hold on
    plot(tv, runoff(:,k), 'r', 'linewidth', lw)
    plot(tv, truth.total_runoff(k,:), 'k', 'linewidth', lw)
    title(['Cell ' num2str(k) ' NSE = ' num2str(nse(k))])
end

% L,T matter. 
% T matters more than L
% Larger L,T are better, but only up to a point
% There are optimal values

figure
for kk=1:n
    subplot(5,4,kk)
    plot(truth.total_runoff(kk,:), runoff(:,kk), 'k.')
    xlabel('Truth')
    ylabel('Predicted')
    title(['Cell ' num2str(kk) ' NSE = ' num2str(nse(kk))])
    hold on
    plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], '--r')
end
