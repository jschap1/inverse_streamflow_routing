% Compare additive Gaussian error vs. multiplicative lognormal error
% 
% Do the Y20 assumptions make A/G error look similar to L/M error?
% 11/5/2023 JRS

clear, clc, close all
% cd /Volumes/HD3/ISR/inverse_streamflow_routing
cd /home/jschap/Documents/ISR/inverse_streamflow_routing/
addpath(genpath('./src/'))

fs=18;lw=2;
load('./allegheny_data/setup/setup-swot-gage.mat');
aa=load('./allegheny_data/setup/setup-2-gage.mat');
basin.distmat = aa.basin.distmat;
load('./allegheny_data/setup/alleg_nldas_nanfilled.mat')
load('./allegheny_data/setup/alleg_tmpa_nanfilled.mat')

n = size(HH,2);
tv = datetime(2009,1,1):datetime(2009,12,31);

%% Load AG errors

runoff_errors = struct();
% this is the runoff prior... need the errors
A = load('./allegheny_data/errors/m0a1L40T5/Ens_errors_AG.mat');
[n, nt, M] = size(A.runoff_errors);
runoff_errors.AG = A.runoff_errors;
clearvars A

%% Load LM errors

B = load('./allegheny_data/errors/m1a1L40T5_LM/simulated_runoff_errors_2000.mat');
runoff_errors.LM = reshape(B.alldata, 20*21, nt, M);
clearvars B
basinmask = basin.mask; % no nan values
basinmask(isnan(basin.mask)) = 0;
basinmask = logical(basinmask);
runoff_errors.LM = runoff_errors.LM(basinmask, :, :); 

%% Compare AG and LM errors

mean(runoff_errors.AG(:))
mean(runoff_errors.LM(:))
std(runoff_errors.AG(:))
std(runoff_errors.LM(:))

%% Conclusion: AG and LM errors have quite different distributions

k = 17; % cell
t = 12; % time
figure
subplot(2,1,1)
histogram(runoff_errors.AG(k,t,:))
title('Add Gaussian errors')
subplot(2,1,2)
histogram(runoff_errors.LM(k,t,:))
title('Mult lognormal errors')

%% AG prior



%% LM prior


