function y = medfilt1_hi(x, Wn, dim)
    n = 2/Wn;
    
    % find closest odd number to n
    if odd(floor(n))
        n = floor(n);
    else
        n = ceil(n);
    end
    if ~odd(n)
        n = n+1;
    end
    
    if nargin < 3
        dim = find(size(x)>1, 1);
    end
    
    y_lo = fastmedfilt1(x, n, dim);    
%     y_lo = medfilt1(varargin{:});
    y = x-y_lo;

end