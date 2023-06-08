% Kalman gain calculation
%
% H = measurement model operator
% P = state error covariance matrix
% R = measurement error covariance matrix
% no_cond = flag to turn off check for ill-conditioned system (uses pinv)
%
% Options for pseudoinverse if the system is ill-conditioned
% Accounts for missing data
% Assumes H, P, R are singles, but could be adapted for sparse matrices

function K = kalman_gain_y20(H, P, R, missing_rows, w)

disp(['Calculating Kalman gain for window ' num2str(w)])

if sum(missing_rows) == 0
%     disp(['Using pinv for window ' num2str(w)])
    K = (P*H')*pinv(H*P*H' + R);
else
    H1 = H;
    H1(missing_rows,:) = 0;
    K = (P*H1')*pinv(H1*P*H1' + R);
%     K = (P*H1')/(H1*P*H1' + R);
    K(isnan(K)) = 0;
end

return