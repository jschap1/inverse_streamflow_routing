% Calculate the standard deviation matrix
%
% This matrix, multiplied element-wise by the correlation matrix, gives us
% the covariance matrix
%
% v = standard deviation of the multiplicative error, epsilon

function phi = calc_phi(x, s, alpha)

[nt,n1] = size(x);

if n1>1
    xflat = reshape(x', n1*(s+1),1);
    phi = alpha^2*xflat*xflat';
else
    phi = alpha^2*x*x';
end

return