% Exceedence ratio
%
% Calculates exceedence ratio for an ensemble prediction according to 
% Moradkhani et al. 2006 paper
%
% truth = true values of variable (n by nt)
% pred = predicted values of variable (n by nt by M, ensemble)
% n = exceedence level (95%, for example)

function er = myER(truth, pred, n)

% we assume a two tailed distribution

upp = prctile(pred,50+n/2,3);
low = prctile(pred,50-n/2,3);

sum1 = sum(truth>upp,2);
sum2 = sum(truth<low,2);

N = sum1+sum2;
T = size(truth,2);
er = N/T;

return