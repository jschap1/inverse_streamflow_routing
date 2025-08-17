% Shaded runoff plot (lognormal)
%
% 8/26/2023 JRS
% Plots runoff central tendency and spread assuming runoff is a lognormally
% distributed random variable

% shadedErrorBar(1:365,m,err,'lineprops','-r','patchSaturation',0.33)

function [m, err] = shaded_runoff_plot_lognorm(runoff_ensemble, cell1)

post_runoff_ensemble_cell = squeeze(runoff_ensemble(cell1,:,:))';
mean_runoff_cell = mean(post_runoff_ensemble_cell);
lpost_runoff_ensemble_cell = log(post_runoff_ensemble_cell);
lmean_cell = mean(lpost_runoff_ensemble_cell);
lsigma_cell = std(lpost_runoff_ensemble_cell);
err1 = exp(lmean_cell - 2*lsigma_cell);
err2 = exp(lmean_cell + 2*lsigma_cell);
err = [err2;err1];
m = mean_runoff_cell;

return