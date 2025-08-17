% Propagates runoff error to discharge error

% g = gauge number
% Uses uncertainty propagation formula (Equation 10 from David et al. 2019)
% Must be done with a monthly timestep so that all flow gets to the outlet
% in one timestep and we can write the network matrix without lag terms.

function stdprop = propagate_runoff_error(stdest, HH, g)

% Create the network matrix
[HH, travel_time] = state_model(basin_mask, basin_mask_gauge, flow_length_gauge, ...
    flow_vel, tstep);

% Network matrix is H, but for m=n (a gauge at every grid cell)
% Also, it is assumed that all runoff gets to outlet in one time step
% So David et al., 2019 worked with monthly errors


nt = size(stdest,1);
[m,n,kp1] = size(HH); % k plus 1

sumvar = zeros(nt,1);
for ki=1:kp1
    ind = logical(HH(g,:,ki));
    sumvar = sumvar + sum(stdest(:,ind).^2,2);
end
stdprop = sqrt(sumvar);

return