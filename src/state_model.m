% calculate all the linear operatores for the state routing model
% F maps runoff to streamflow over all points
% H maps runoff to streamflow over gauges only

function [HH, travel_time] = state_model(basin_mask, basin_mask_gauge, flow_length_gauge, ...
    flow_vel, tstep)

[nrows ncols] = size(basin_mask);
ncells = nansum(nansum(basin_mask));

if max(flow_length_gauge)==0
    ksteps = ceil(max(max(max(1)/flow_vel/tstep)));
else
    ksteps = ceil(max(max(max(flow_length_gauge)/flow_vel/tstep)));
end

ngauges = size(flow_length_gauge, 3);

HH = zeros(ngauges, ncells, ksteps);

% create the mapping for compacting cells
% tmp_mask_flat = basin_mask(:);
% tmp_mask_flat(isnan(tmp_mask_flat)) = 0;
% cmap = cumsum(tmp_mask_flat);
% gauge_loc = cmap(gauge_i+(gauge_j-1)*nrows);

travel_time = cell(ngauges,1);
for g=1:ngauges
%
%     % Find out which cell(s) are not contributing and handle them
%     contributing_mask = squeeze(flow_length_gauge(:, :, 1));
%     contributing_mask(~isnan(contributing_mask)) = 1;
%     contributing_mask(isnan(contributing_mask)) = 0;
%     bmask = basin_mask;
%     bmask(isnan(bmask)) = 0;
%     figure,imagesc(bmask - contributing_mask)
%     figure,imagesc(basin_mask)

    flow_length_gauge_compact = reshape(squeeze(flow_length_gauge(:, :, g)), nrows*ncols, 1);
    flow_length_gauge_compact(isnan(flow_length_gauge_compact)) = [];

    travel_time{g} = flow_length_gauge_compact/flow_vel/tstep;
    travel_steps = floor(travel_time{g}) + 1; %Y

    % Try one of these:
%     travel_steps = round(travel_time{g}) + 1;%M
%     travel_steps = floor(travel_time{g}) + 1; %Y

    basin_mask_gauge_compact = reshape(squeeze(basin_mask_gauge(:, :, g)), nrows*ncols, 1);
    basin_mask_gauge_compact(isnan(basin_mask_gauge_compact)) = [];

    for s=1:ksteps
        HH(g, travel_steps==s, s) = 1;
        HH(g, basin_mask_gauge_compact~=1, s) = 0;
    end

end

return
