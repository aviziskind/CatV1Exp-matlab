function y = wgt_nanstd(x,w, x_cent)
    if (nargin < 2) || isempty(w)
        y = nanstderr(x);
        return;
    end        
    i = ~isnan(x);
    x=x(i);
    w=w(i);
    x = x(:); w = w(:);            
    w = w/sum(w);
    
    if nargin < 3
        x_cent = sum(x(:) .* w(:));
    end
    
    x_diff = x - x_cent;
    
    var_x = sum( w .* (x_diff.^2) );
    y = sqrt(var_x);
end

% formula for weighted variance:
%     w = w/sum(w);
%     wgt_mn = sum(w.*x);
%     x_d = x - wgt_mn;
%     var = sum(w .* x_d.^2)
