% Calculated augmented model operators
%
% 9/22/2020 JRS
% Ming's approach
% Sparse matrix option

function [Hprime, L] = calc_augmented_model(H, sm)

[ngage, ncells, lag_steps] = size(H);

if issparse(H)
    Hprime = sparse(sm*ngage, sm*ncells);
else
    Hprime = zeros(sm*ngage, sm*ncells);
end

% Reshaping H
H_flat = reshape(H, ngage, ncells*lag_steps);

% Create the upper diagonal Hprime matrix (code copied from Ming)
% (have tested this; it works)
for s=1:sm
    if s+lag_steps <= sm
        Hprime(ngage*(s-1)+1:ngage*s, ncells*(s-1)+1:ncells*(s-1+lag_steps)) = H_flat;
    else
        Hprime(ngage*(s-1)+1:ngage*s, ncells*(s-1)+1:ncells*sm) = H_flat(:, 1:ncells*(sm-s+1));
    end
end

% Create the lower diagonal L matrix

if issparse(H)
    L = sparse(sm*ngage, sm*ncells);
else
    L = zeros(sm*ngage, sm*ncells);
end

empty_rows = sm + 1 - lag_steps;
empty_cols = empty_rows;
m = ngage;
n = ncells;

% k = lag_steps - 1; % using Pan and Wood (2013) notation

for s_col = (empty_cols+1):sm
    
    d2 = ((s_col-1)*n + 1):s_col*n; % column index (block matrix)
    
    % step through columns of the block matrix
    for s_row=(empty_rows+1):sm
        
        d1 = ((s_row-1)*m+1):s_row*m; % row index (block matrix)
        
        if s_row >= s_col
            ind = lag_steps - (s_row - s_col);
            if ind > 0
                L(d1,d2) = H(:,:,ind);
            end
        end

    end
end


return

%% More scrap

% % Populate H matrix by blocks
% for r = 1:(s+1)
%     for c=1:(s+1)
%         
%         r_ind = ((r-1)*n+1):r*n;
%         c_ind = ((c-1)*m+1):c*m;
%         
%         if c<(k+1)
%             H(r_ind,c_ind) = HH(:,:,c);
%         end
%         
%         
%         
%     end
% end
% 
% 
% blkdiag(HH(:,:,1), 3)
% 
% blkdiag(HH(:,:,1), HH(:,:,1), HH(:,:,1), HH(:,:,1), HH(:,:,1), HH(:,:,1))

%% Scrap
  
%         % populate the diagonal of L
%         if s_col == s_row
%             L(d1, d2) = H(:,:,4);
%         end
%         
%         % one diag to the left (for all rows below the first non-empty row)
%         if s_row - s_col == 1 && s_row 
%             L(d1,d2) = H(:,:,3);
%         end
% %         
% %         % two diags to the left
% %         if s_row - s_col == 2
% %             L(d1,d2) = H(:,:,2);
% %         end
%     
%     end
%     
% end
% 
% 
% 
% 
% 
% 
% % Create the lower diagonal L matrix
% L = zeros(sm*ngage, sm*ncells);
% 
% % do the first (m by n) block
% m = ngage;
% n = ncells;
% L(1:m, 1:n) = zeros(m,n); % 1:2, 1:3
% 
% % next block
% L((m+1):2*m, (n+1):2*n) = zeros(m,n); % 3:4, 4:6
% 
% % continue until end of L
% L((2*m+1):3*m, (2*n+1):3*n) = zeros(m,n); % 5:6, 7:9
% 
% % Do the first row of blocks (all zeros)
% L = zeros(sm*ngage, sm*ncells);
% for s=1:sm
%     d1 = 1:m;
%     d2 = ((s-1)*n + 1):s*n;
%     L(d1,d2) = ones(m,n);
% end
% 
% %     d1 = ((s-1)*m + 1):s*m;
% %     d2 = 1:n;
% 
% % Do the second row of blocks
% for s=1:sm
%     d1 = (m+1):2*m;
%     d2 = ((s-1)*n + 1):s*n;
%     L(d1,d2) = ones(m,n);
% end
% 
% % Third row
% for s=1:sm
%     d1 = (2*m+1):3*m;
%     d2 = ((s-1)*n + 1):s*n;
%     L(d1,d2) = ones(m,n);
% end
% 
% % Fourth row
% for s=1:sm
%     d1 = (3*m+1):4*m;
%     d2 = ((s-1)*n + 1):s*n;
%     L(d1,d2) = ones(m,n);
% end
% 
% % Fifth (last) row
% for s=1:sm
%     d1 = (4*m+1):5*m;
%     d2 = ((s-1)*n + 1):s*n;
%     L(d1,d2) = ones(m,n);
% end
% 
% H = ones(ngage, ncells, lag_steps);
% L = zeros(sm*ngage, sm*ncells);
% 
% empty_rows = sm + 1 - lag_steps;
% % figure,imagesc(L)
% % once we have reached the (empty_rows + 1) row, we begin populating L
% 
% % Third row (the first non-empty row of the block matrix)
% for s=1:sm
%     
%     d2 = ((s-1)*n + 1):s*n;
%     
%     % step through columns of the block matrix
%     for s_col = 1:sm
%         d1 = ((s_col-1)*m+1):s_col*m;
%         
%     d1 = (2*m+1):3*m;
%     
%     
%     if s==2 % second block in the row
%         L(d1,d2) = H(:,:,4); % H for lag3 (the last lag)
%     end
%     
% end
% 
% % Fourth row
% for s=1:sm
%     d1 = (3*m+1):4*m;
%     d2 = ((s-1)*n + 1):s*n;
%     
%     if s==2 % second block
%         L(d1,d2) = H(:,:,3);
%     elseif s==3 % third block
%         L(d1,d2) = H(:,:,4);
%     end
%     
% end
% 
% % Fifth (last) row
% for s=1:sm
%     d1 = (4*m+1):5*m;
%     d2 = ((s-1)*n + 1):s*n;
%     
%     if s==2 % second block
%         L(d1,d2) = H(:,:,2);
%     elseif s==3 % third block
%         L(d1,d2) = H(:,:,3);
%     elseif s==4 % fourth block
%         L(d1,d2) = H(:,:,4);
%     end
%     
% end
% 
% 
% % empty for the first (sm - 1 + lag_steps) rows
% 
% 
% 
% for s=1:sm
%     if s+lag_steps <= sm
%         L(ngage*(s-1)+1:ngage*s, ncells*(s-1)+1:ncells*(s-1+lag_steps)) = H;
%     else
%         L(ngage*(s-1)+1:ngage*s, ncells*(s-1)+1:ncells*sm) = H(:, 1:ncells*(sm-s+1));
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% % L is empty for the first ngage*(sm - lag_steps + 1) rows 
% % and ncells(sm - lag_steps + 1) columns
% i1 = ngage*(sm - lag_steps + 1);
% j1 = ncells(sm - lag_steps + 1);
% 
% % after that, the ngage*(lag_steps - 1) by ncells*(lag_steps - 1) block is
% % filled with H1, H2, ... H_(lag_steps - 1) in the lower diagonal.
% 
% [ni,nj] = size(L);
% for i=1:ni
%     for j=1:nj
%         if i <= i1 || j <= j1
%             continue
%         else
%             L(i,j) = 1;
%         end
%     end
% end
% 
% % augmented model operator Hprime
% H0 = H(:,:,1);
% H1 = H(:,:,2);
% H2 = H(:,:,3);
% H3 = H(:,:,4);
% Z = zeros(size(H0));
% 
% L1 = [Z, Z, Z, Z;
%     Z, H3, Z, Z;
%     Z, H2, H3, Z;
%     Z, H1, H2, H3];
% 
% Hprime1 = [H0, H1, H2, H3;
%     Z, H0, H1, H2;
%     Z, Z, H0, H1;
%     Z, Z, Z, H0];
% 
% Hprime(1,1:20) = H0;
% Hprime(2,21:40) = H1;
% Hprime(3,41:60) = H2;
% Hprime(4,61:80) = H3;
% 
% 
% L(2,21:40) = H3;
% L(3,41:60) = H2;
% L(3,61:80) = H3;
% L(4,21:40) = H1;
% L(4,41:60) = H2;
% L(4,61:80) = H3;
% 
% 
% 
% 
% return