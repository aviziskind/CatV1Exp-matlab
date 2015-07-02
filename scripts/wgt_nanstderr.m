function y = wgt_nanstderr(x,w)
    if (nargin < 2) || isempty(w)
        y = nanstderr(x);
        return;
    end        
    i = ~isnan(x);
    x=x(i);
    w=w(i);
    x = x(:); w = w(:);    
    w = w/sum(w);

    wgt_x = x(:) .* w(:);  wgt_nanstderr    
    x_diff = x - sum(wgt_x);
    
    var_x = sum( w .* (x_diff.^2) );
    n = length(x);
    y = sqrt(var_x / n);
end