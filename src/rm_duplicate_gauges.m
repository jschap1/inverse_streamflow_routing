function [gage_lon, gage_lat, gage, bb] = rm_duplicate_gauges(gage_lon, gage_lat, gage)

gage_sums = nansum(gage,2);

% test case
% gage_lon = [4;5;6;6];
% gage_lat = [1;1;1;1];
% gage_sums = [4,5,6,7]';

% Keep downstream-most gauges
[aa,bb,cc] = unique([gage_lon, gage_lat], 'rows');
ngauges = length(bb); % number of unique gauges
gages_to_keep = zeros(ngauges,1);
for kk=1:ngauges
    ind = find(cc==kk); % find the gauge that is most downstream within a cell
    [~, maxind] = max(gage_sums(ind));
    gages_to_keep(kk) = ind(maxind);
end

% [gage_lon(cc(ind)), gage_lat(cc(ind))]
% gage_sum(cc(ind))

gage = gage(gages_to_keep,:);
gage_lon = aa(:,1);
gage_lat = aa(:,2);

return