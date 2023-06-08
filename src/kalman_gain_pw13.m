% Kalman gain calculation
%
% H = measurement model operator
% P = state error covariance matrix
% R = measurement error covariance matrix
% no_cond = flag to turn off check for ill-conditioned system (uses pinv)
%
% Options for pseudoinverse if the system is ill-conditioned
% Accounts for missing data

function K = kalman_gain_pw13(H, P, R, missing_rows, w)

if sum(missing_rows) == 0
    disp(['Using pinv for window ' num2str(w)])
    % using sparse matrix: 257 seconds per window for m=240, s=2*(k+1), n = 3681
    % using sparse matrix: 52 seconds per window for m=240, s=k+1, n = 3681
    K = (P*H')/(H*P*H' + R);
%     K = (P*H')*pinv(H*P*H' + R);
else
    H1 = H;
    H1(missing_rows,:) = 0;
    K = (P*H1')/(H1*P*H1' + R);
    K(isnan(K)) = 0;
end

return