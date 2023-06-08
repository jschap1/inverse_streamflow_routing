% Inverse routing
% 
% 4/17/2023
% EnKF version of ISR, with domain localization
% Functions as a wrapper for ISR_EnKF, calling it for each subdomain
%
% INPUTS
% runoff_prior = prior runoff ensemble (n, nt, M)
% HH = measurement model operator
% discharge = gage measurements
% s = length of smoothing window
%
% OUTPUTS
% runoff = posterior runoff ensemble
% w1 = start times for each window (for the last subbasin)
% w2 = end times for each window (for the last subbasin)
% maxv = the largest state vector used in the domain-ISR

function [runoff, w1, w2, maxv, maxv_sub] = ISR_domain_wrapper(runoff_prior, HH, gage, s, basin, Cv, opt, tv)

maxv_sub = 0; % initialize

% Sort gauges from upstream to downstream (by drainage area)
[nt,m] = size(gage);
M = size(runoff_prior,3);
[gage_area_sorted,ind_inc_DA] = sort(basin.gage_area); % gauges, sorted from upstream to downstream
n = size(HH,2);

% runoff = zeros(n,nt,M);
opt.s_total = s;

sm = 1e-3; % a small positive value to replace spurious negative values in the adjusted discharge

% For each gauge, make a list of upstream gauges
% identify the cell the upstream gauges are in
us_gages = cell(m,1);
us_gages_tt = cell(m,1);
cell_inds = (1:n)';
gage_inds = sub2ind(size(basin.mask), basin.gage_i, basin.gage_j);
basin.mask(isnan(basin.mask)) = 0;
basin.mask = logical(basin.mask);

% Calculate travel time to each gage (for later use)
[HH, travel_time] = state_model(basin.mask, basin.mask_gage, basin.flow_length_gage, ...
    basin.flow_vel, basin.timestep);
opt.HH_total = HH;
opt.domISR = 1;

cells_run = NaN(size(basin.mask));
runoff = runoff_prior; % initialize

for mm=1:m
    
    if mm==93
        1
    end
    
    current_gage = ind_inc_DA(mm);
    
    ttmap = flipud(make_map(basin, travel_time{current_gage}));
     % travel time from this gauge to the outlet
    us_gages_tt{current_gage} = ttmap(gage_inds); % travel time to outlet of current subbasin from each upstream gauge
    % value is zero if the gauge is not connected
    
    us_gages{current_gage} = find(us_gages_tt{current_gage}); 
    
%     us_gages_tt{current_gage} = ceil(us_gages_tt{current_gage}); % round up to the next whole time step
    
    us_gages_tt{current_gage} = round(us_gages_tt{current_gage}); % not used?
    
%     us_gages_tt{current_gage} = round(us_gages_tt{current_gage});
    
    
    % Get the grid cells that are part of the current sub-basin
    current_mask = basin.mask_gage(:,:,current_gage);
    current_mask_all = current_mask;
    
    % remove any grid cells that have already been accounted for
    current_mask(cells_run==1) = 0;

    current_mask_lin = current_mask(basin.mask);
    current_mask_all_lin = current_mask_all(basin.mask);
    ind_current_basin = find(current_mask_lin);
    ind_current_basin_all = find(current_mask_all_lin);
    cells_run_lin = cells_run(basin.mask);
    ind_cells_run = find(cells_run_lin==1);
    
    % Calculate local state model

    % crop to the smallest possible rectangle
    [keeprows, keepcols] = find(current_mask==1);
    keeprows = unique(keeprows);
    keepcols = unique(keepcols);
    basinmask = current_mask(keeprows,keepcols);
    basinmask(basinmask == 0) = NaN;
    flow_length_gauge = basin.flow_length_gage(keeprows,keepcols,current_gage);
    flow_length_gauge(isnan(basinmask)) = NaN;
    basin_mask_gauge = basin.mask_gage(keeprows,keepcols,current_gage);
    basin_mask_gauge(isnan(basinmask)) = NaN;
    
    [HH_sub, ~] = state_model(basinmask, basin_mask_gauge, ...
        flow_length_gauge, basin.flow_vel, basin.timestep);
    
    k_sub = size(HH_sub,3) - 1;
    
%     s_sub = 2*(k_sub+1); % bigger K
    s_sub = k_sub+1;
    
    Q = gage(:,current_gage);
%     Q_orig = Q;
       
%     figure
%     plot(Q, 'linewidth', 2)
%     title(['Gage ' num2str(current_gage)])
%     ylabel('discharge (mm/hr)')
%     xlabel('time')
    
    opt.current_gage = current_gage;
    opt.ind_current_basin = ind_current_basin;
    [runoff, w1, w2] = ISR_domain(runoff, HH_sub, Q, s_sub, basin, Cv, opt);    
    
    % check results
%     figure
%     lw=2;
%     plot(runoff_prior(ind_current_basin, :, 1), 'blue', 'linewidth', lw)
%     hold on
%     plot(runoff(ind_current_basin, :, 1), 'red', 'linewidth', lw)
%     plot(basin.true_runoff(ind_current_basin,:), 'k', 'linewidth', lw)
%     legend('prior','posterior', 'true')
    
    cells_run(current_mask==1) = 1;
    
    if mm==50
        1
    end
    
    n_sub = length(ind_current_basin);
    if n_sub*(s_sub+1)>maxv_sub
        maxv_sub = n_sub*(s_sub+1);
    end
    
    disp(['Finished sub-basin ' num2str(mm) ' of ' num2str(m)])
    
end

maxv = n*(s+1);

disp(['Domain localization reduced the length of the state vector by ' num2str(100-100*maxv_sub/maxv,3) '%'])
% figure,imagesc(mean(runoff,3)), colorbar, title('est runoff'), caxis([0,7])
% figure,imagesc(basin.true_runoff), colorbar, title('true runoff'), caxis([0,7])

return