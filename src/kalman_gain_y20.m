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

function K = kalman_gain_y20(H, P, R, missing_rows, w, varargin)

if nargin>5
    sparseflag = varargin{1};
    disp('Using sparse matrices for kalman gain')
else
    sparseflag = 0;
    disp('Using singles for kalman gain')
end

switch sparseflag
    
    case 1
        
        % Make any single matrices sparse
        if ~issparse(P)
            P = sparse(P);
        end
        if ~issparse(H)
            H = sparse(H);
        end
        if ~issparse(R)
            R = sparse(R);
        end
        
    case 0
        
        % Make any sparse matrices full
        if issparse(P)
            P = single(full(P));
        end
        if issparse(H)
            H = single(full(H));
        end
        if issparse(R)
            R = single(full(R));
        end
        
end

% Cut to size for debugging
% m = 240;
% n = 3681;
% H = H(1:m,1:n);
% P = P(1:n,1:n);
% R = R(1:m,1:m);
% missing_rows = missing_rows(1:m);

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

% innovation = NaN(m,1);
% x_update = K*innovation;

return