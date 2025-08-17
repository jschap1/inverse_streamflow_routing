% Does the reverse operation of fliparoo

function x = unfliparoo(xw, n, s)
    
    x = flipud(reshape(xw, n, s+1)');
    
return

