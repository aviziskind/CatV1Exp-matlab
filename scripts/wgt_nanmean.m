function m = wgt_nanmean(x,w)
    if (nargin < 2) || isempty(w)
        m = nanmean(x);
        return;
    end    

    idx = ~isnan(x);
    x=x(idx);
    w=w(idx);
    
    if isempty(x)
        m = nan;        
    else    
        w = w/sum(w);
        m = sum(x(:) .* w(:));
    end
end