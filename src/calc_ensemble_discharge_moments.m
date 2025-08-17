% Calculates ensemble discharge and moments

function [meanQ, stdQ, Qens,medianQ,qlow, qupp, Q_IQR] = calc_ensemble_discharge_moments(runoff, HH)

[n, nt, M] = size(runoff);
[ngage, n, kp1] = size(HH);
Qens = zeros(nt, ngage, M);
for mm=1:M
    Qens(:,:,mm) = state_model_dumb(runoff(:,:,mm)', HH);
end
meanQ = mean(Qens,3);
stdQ = std(Qens,[],3);

medianQ = median(Qens,3);
qlow = prctile(Qens, 25 ,3);
qupp = prctile(Qens, 75 ,3);

Q_IQR = iqr(Qens,3); % goes with median (acts like stddev)

% figure,plot(stdQ(:,1))
% hold on, plot(r(:,1))
% legend('stddev','IQR')
% % 129 to 206 IQR for Q(20,1,:) prior
% % 174 to 218 for posterior -> tighter bounds. good.
return