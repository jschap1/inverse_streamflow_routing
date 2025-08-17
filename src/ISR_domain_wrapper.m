% Inverse routing
% 
% 4/17/2023
% EnKF version of ISR, with domain localization
% Functions as a wrapper for ISR_EnKF, calling it for each subdomain
% Can also be used with PW13 or Y20 ISR. Just need to enable the flag.
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

Y20_flag = 0; % flag for using Y20 ISR

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
opt.gages = gage;

cells_run = NaN(size(basin.mask));
runoff = runoff_prior; % initialize

% list of gauges that have been upstream gauges
% these are removed from the upstream gauge list to avoid double counting
% when doing the gauge correction for upstream subdomains
usglist = zeros(0,1);

for mm=1:m
    
    % breakpoint to see what is going on within a particular subbasin
    if mm==15
        1
    end
    
    current_gage = ind_inc_DA(mm);
    
    ttmap = flipud(make_map(basin, travel_time{current_gage}));
    ttmap = make_map(basin, travel_time{current_gage});
    
    opt.plot=0;
    if opt.plot
        figure
        plotraster(basin.lonv, basin.latv, flipud(ttmap), 'tt')
        hold on
        plot(basin.gage_lon, basin.gage_lat, 'r.', 'markersize', 20)
    end
    
     % travel time from this gauge to the outlet
    us_gages_tt{current_gage} = ttmap(gage_inds); 
    % travel time to outlet of current subbasin from each upstream gauge
    % value is zero if the gauge is not connected
    
    
    
%     us_gages_tt{current_gage} = ceil(us_gages_tt{current_gage}); % round up to the next whole time step
    
    % we need this in order to subtract upstream discharge
    us_gages{current_gage} = find(us_gages_tt{current_gage}); 
    us_gages_tt{current_gage} = round(us_gages_tt{current_gage});
    % we need only the immediate upstream gauges, not additional gauges
    % beyond that
    
    
    
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
    
    % get distance matrix for current subbasin
    subbasin = struct();
    subbasin.grid_lon = basin.grid_lon(keeprows, keepcols);
    subbasin.grid_lat = basin.grid_lat(keeprows, keepcols);
    subbasin.lon = subbasin.grid_lon(~isnan(basinmask));
    subbasin.lat = subbasin.grid_lat(~isnan(basinmask));
    subbasin.distmat = calc_distance_matrix(subbasin.lon, subbasin.lat);
    
    [HH_sub, ~] = state_model(basinmask, basin_mask_gauge, ...
        flow_length_gauge, basin.flow_vel, basin.timestep);
    
    k_sub = size(HH_sub,3) - 1;
    
    s_sub = 2*(k_sub+1); % bigger K
%     s_sub = k_sub+1;
    
    % Discharge at current gauge
    Q = gage(:,current_gage);
    
    % Discharge at current gauge minus (posterior) time-lagged discharge at
    % upstream gauges
%     nupstream_gages = length(us_gages{current_gage});
    current_upstream_gages = us_gages{current_gage};
    
    % remove any upstream gauges that have already been counted
    current_upstream_gages = current_upstream_gages(~ismember(current_upstream_gages,usglist));
%     if nupstream_gages>0 % if there are upstream gauges
%         meas_ind = ~isnan(Q);
%         lagtime = us_gages_tt{current_gage};
%         Qcorrected = Q;
%         for gind = 1:nupstream_gages
%             lagtime(current_upstream_gages(gind))
%             Qcorrected(meas_ind) = Q(meas_ind) - gage_updated(meas_ind,current_upstream_gages(gind));
%         end
%     end
%     Problem with subtracting upstream Q is that with meas errors,
%     upstream Q may be larger than ds Q. Worse if Q is already updated.
    
    
%     Q_orig = Q;
       
%     figure
%     plot(Q, 'linewidth', 2)
%     title(['Gage ' num2str(current_gage)])
%     ylabel('discharge (mm/hr)')
%     xlabel('time')
    
    opt.current_gage = current_gage;
    opt.current_upstream_gages = current_upstream_gages;
    opt.us_gages_lag = us_gages_tt{current_gage};
    opt.ind_current_basin = ind_current_basin;
    if Y20_flag
        % Y20 ISR
        opt.use_domain_loc = 1;
        % need to pass in distance matrix for subdomain
        basin.distmat = subbasin.distmat;
        runoff = ISR_Y20_domain(runoff, HH_sub, Q, s_sub, basin, ...
            opt, opt.L, opt.T, opt.rho_thres);
    else
        % Ensemble ISR
        [runoff, w1, w2] = ISR_domain(runoff, HH_sub, Q, s_sub, basin, Cv, opt);    
        gage_updated = state_model_dumb(mean(runoff,3)', HH); % using ensemble mean
    end
    
    usglist = [usglist; us_gages{current_gage}]; % add upstream gauges to usglist
    usglist = unique(usglist);
    
    % check results
    
    opt.plot=0;
    if opt.plot
        % does posterior discharge match the gauge?
        ypost = state_model_dumb(mean(runoff,3)', HH);
        figure
        plot(ypost(:,current_gage))
        hold on
        plot(gage(:,current_gage), 'k.', 'markersize', 20)
    end
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