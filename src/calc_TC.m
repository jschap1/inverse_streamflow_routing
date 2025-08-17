% Calculate temporal correlation matrix

function TC = calc_TC(k2, rho_thres, n, s, RAM)

% disp('Calculating temporal correlation matrix')
% note: the minimum size for TC can be quite large even for 
% small T bc of the block diagonal of ones

% determine the proportion of nonzero elements in the TC matrix
maxlag = -log(rho_thres)/k2; 
% floor(maxlag)*2+1; % number of nonzero block bands in the TC matrix
nnzblocks = ((s+1)-2)*(floor(maxlag)*2+1) + 2*(floor(maxlag)+1);
nblocks = (s+1)^2;
percent_nonzero = nnzblocks/nblocks;

% once the sparse matrix has more nonzero elements than this (by percent)
% a single precision matrix will use less memory
percent_nonzero_thres = (n*(s+1)*n*(s+1)-2*n*(s+1)-2)/(4*n*(s+1)*n*(s+1));

% memory required for the sparse matrix
nnz = percent_nonzero*(n*(s+1))^2;
D = max(nnz,1)*16+8*(n*(s+1)+1);
D = D/1e9; % memory requirement (GB)
% disp(['Sparse matrix will require ' num2str(D) ' GB of memory'])

% memory required for the single matrix
S = 4*((n*(s+1))^2)/1e9;
% disp(['Single-precision matrix will require ' num2str(S) ' GB of memory'])

if percent_nonzero < percent_nonzero_thres
    % use sparse
%     disp('Using sparse matrix for TC')
    TC = sparse(n*(s+1),n*(s+1));
else
    % use single
%     disp('There are many non-zero values in TC')
%     disp('Using single precision matrix for TC')
    TC = zeros(n*(s+1),n*(s+1), 'single');
end

% n_thres = 500; % maximum size beyond which we use sparse matrices

% Calculate time lag matrix
% ll = cell(s+1,1);
% for t=1:(s+1)
%     ll{t} = (t-1)*ones(n,n);
% end

% Trying a more efficient method
ll = cell(s+1,1);
for t=1:(s+1)
    ll{t} = (t-1);
end

% % Assemble augmented time-lag matrix
% if n<=n_thres
%     L1 = zeros(n*(s+1),n*(s+1));
% %     L1 = zeros(n*(s+1),n*(s+1), 'int8');
% else
%     L1 = sparse(n*(s+1),n*(s+1));
% %     L1 = zeros(n*(s+1),n*(s+1), 'int8');
% end

% % Assemble augmented time-lag matrix
% for r=1:(s+1)
%     for c=1:(s+1)
%         L1((r-1)*n+1:r*n, (c-1)*n+1:c*n) = ll{abs(r-c)+1};
%     end
%     disp(r);
% end

% Eigenvalue decomposition (L1 is too big to do eigenvalue decomposition)
% [V,D] = eig(L1);
% L2 = V*D*inv(V); % should match L1

% rho_thres = exp(-3); % affects how large the sparse matrix is

 % Calculate TC without storing L1
 for r=1:(s+1)
    for c=1:(s+1)
        vals = exp(-k2*ll{abs(r-c)+1});
        vals(vals<rho_thres) = 0;
        % note: sparse is only supported for double and logical
        
        TC((r-1)*n+1:r*n, (c-1)*n+1:c*n) = vals;
        % note: indexing into sparse matrix is bad practice
        % see calc_SC for a better way
    end
%     disp(r)
 end
    
% Calculate temporal correlation matrix
% TC = exp(-k2*L1); % cannot compute exp of an integer 
% also, the result, will be small

% TC(TC<rho_thres) = 0;

% figure
% imagesc(TC)

return