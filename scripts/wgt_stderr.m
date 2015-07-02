function m = wgt_stderr(x,w)
    if (nargin < 2) || isempty(w)
        w = ones(size(x));
    end        
    i = ~isnan(x);
    x=x(i);
    w=w(i);
    m = stderr(x .* w);
end