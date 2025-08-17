% Extracts runoff for a particular subdomain in a basin
%
%

function sub1runoff = get_runoff_for_cells_in_subdomain(runoff, subdomain1_index, basin)

[nr, nc, ~] = size(basin.mask_gage);
[nt,n] = size(runoff);
n1 = length(subdomain1_index);
tmpa_runoff_prior_map = zeros(nr,nc,nt);
sub1runoff = zeros(n1,nt);
for t=1:nt
    tmpa_runoff_prior_map(:,:,t) = make_map(basin, runoff(t,:));
    tmp = tmpa_runoff_prior_map(:,:,t);
    sub1runoff(:,t) = tmp(subdomain1_index);
end

return