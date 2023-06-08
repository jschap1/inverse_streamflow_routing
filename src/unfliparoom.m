% Does the reverse operation of fliparoo
% 
% For multiple inputs (ensemble case)

function x = unfliparoom(xw, n, s)
    
    M = size(xw, 2);
    x = zeros(s+1,n,M);
    for i=1:M
        x(:,:,i) = flipud(reshape(xw(:,i), n, s+1)');
    end
    
return

