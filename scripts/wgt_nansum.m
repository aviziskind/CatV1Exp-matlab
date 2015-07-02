function S = wgt_nansum(x,w)
    if (nargin < 2) || isempty(w)
        m = nansum(x);
        return;
    end    

    idx = ~isnan(x);
    x=x(idx);
    w=w(idx);
    w = w/sum(w);
    S = sum(x(:) .* w(:));
    
end