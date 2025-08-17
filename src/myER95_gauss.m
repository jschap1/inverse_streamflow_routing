% Calculates 95% exceedence ratio for Gaussian data

function er95 = myER95_gauss(truth, pred, predsd)

top5 = pred + 1.96*predsd;
bottom5 = pred - 1.96*predsd;
sum1 = sum(truth>top5,2);
sum2 = sum(truth<bottom5,2);
T = size(top5,2);
er95 = (sum1+sum2)/T;

return