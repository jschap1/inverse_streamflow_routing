% Calculate spatial correlation matrix

function SCaug = calc_SC(k1, distmat, rho_thres, n, s, RAM)

% disp('Calculating spatial correlation matrix')

% determine the proportion of nonzero elements in the SC matrix
maxdist= -log(rho_thres)/k1;
nnz = sum(distmat(:)<=maxdist);
% floor(maxlag)*2+1; % number of nonzero bands in the SC matrix
n_entries = (n)^2;
percent_nonzero = nnz/n_entries;

% once the sparse matrix has more nonzero elements than this (by percent)
% a single precision matrix will use less memory
% percent_nonzero_thres = (n*(s+1)*n*(s+1)-2*n*(s+1)-2)/(4*n*(s+1)*n*(s+1));

% memory required for the sparse matrix
nnz = percent_nonzero*(n*(s+1))^2;
D = max(nnz,1)*16+8*(n*(s+1)+1);
D = D/1e9; % memory requirement (GB)
% disp(['Sparse matrix will require up to ' num2str(D) ' GB of memory'])

% memory required for the single matrix
S = 4*((n*(s+1))^2)/1e9;
% disp(['Single-precision matrix will require ' num2str(S) ' GB of memory'])

if D>=RAM && S>=RAM
    
    error('SCaug matrix will not fit in memory')
    
elseif D <= S
    
    % use sparse
%     disp('Using sparse matrix for SC')
        
    SC = exp(-k1*distmat);
    SC(SC<rho_thres) = 0;
    
    % check nonzero fraction
%     sum(SC(:)>0)/numel(SC);
        
    % Find indices of nonzero values
    [i,j] = find(SC~=0);
    lin_ind = logical(SC(:)~=0);
    
%     i_aug = zeros(nnz,1);
%     j_aug = zeros(nnz,1);
%     for ii=1:(s+1)
%         for jj=1:(s+1)
%             i_aug(ii)
%             j_aug(jj)
%         end
%     end
 
%     i = [1;2;2];
    
    % Create i_aug and j_aug, the indices in SCaug 
    % where the nnz nonzero values are found
%     i_aug = ones(nnz*(s+1)^2,1);
    
    % nnz = number of nonzero values in the SC matrix
    nnz = sum(distmat(:)<=maxdist);

    % Create i_aug
    i_aug = zeros(nnz*(s+1)^2,1);
    ind = 1:length(i)*(s+1);
    for ss=0:s
        i_aug(ind) = repmat(i, s+1, 1) + n*repelem(0:s, nnz)';
        ind = ind + (s+1)*nnz; % should it be ss*nnz?
%         disp(ss)
    end % took 37.5 s for ORB
    
    % Create j_aug
    j_aug = zeros(nnz*(s+1)^2,1);
    ind = 1:length(j)*(s+1);
    for ss=0:s
        j_aug(ind) = repmat(j + ss*n, s+1, 1);
        ind = ind + (s+1)*nnz;
%         disp(ss)
    end % took 30.4 s for ORB
    
    % i_aug and j_aug and v_aug require a lot of RAM (10GB each)
    % To store SCaug (22GB) at the same time will be too much
    
    
    
    % create i_aug, j_aug
    % inputs: n, s, i, j
%     j_aug = repmat(j, (s+1)^2, 1)
%     u = zeros(nnz*(s+1)^2,1);
    % s+1 or nnz? 
    % need a better example w/out nnz = s+1 = n
    
%     i_aug = repmat(i, s+1, 1);
    
    
    % Need to add n to the index, every s+1 steps
    % to populate the whole SCaug matrix
%     v = n*[0:s];
%     u = repelem(v,length(i));
%     i_aug = i_aug + u';
%     j_aug = j_aug + u';
    
    % Make nonzero value array
    v_aug = repmat(SC(lin_ind), (s+1)^2, 1);
    SCaug = sparse(i_aug, j_aug, v_aug, n*(s+1), n*(s+1));

    % for some reason, when I calculate SCaug this way, there are 32 times
    % more zeros that there should be...
    % This has to do with v_aug
    % There are fewer than expecetd v_aug values
    % Why is this? 
    % Come back to this tonight/tomorrow (10/7/22)
    
%     try
%         SCaug = spalloc(n*(s+1), n*(s+1), nnz);
%     catch
%         error('SCaug is too large')
%         SCaug = sparse(n*(s+1),n*(s+1));
%     end
    
elseif S < D
    
    % use single
%     disp('There are many non-zero values in SC')
%     disp('Using single precision matrix for SC')
    
    SC = single(exp(-k1*distmat));
    SC(SC<rho_thres) = 0;
    SCaug = zeros(n*(s+1),n*(s+1), 'single');
    
    for kk=1:(s+1)
        for jj=1:(s+1)
            SCaug(((kk-1)*n + 1):kk*n, ((jj-1)*n + 1):jj*n) = SC;
        end
    end

end

return