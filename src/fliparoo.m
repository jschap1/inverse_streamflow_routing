% Rearranges the runoff input from (nt by n) to the runoff vector for the
% current ISR window
%
% unfliparoo reverses this

function xw = fliparoo(x, s)
    
    n = size(x,2); % number of cells
    xw = flipud(x);
    xw = reshape(xw', n*(s+1), 1);
    
return

