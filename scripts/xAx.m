function y = xAx(A,x, m)
    % performs the operation matrix operation x'*A*x, for (columns of) vector(s) x, and matrix A. if
    % a third argument, m, is supplied, the result is (x-m)'*X*(x-m)

    d = size(A,1);
    assert(size(A,2) == d);
    if nargin > 2
        assert(all(size(m) == [d,1]));    
        x = bsxfun(@minus, x, m);
    end
    
    [d_x, N] = size(x);
    if (d_x ~= d) && (N == d);
        x = x';
        [d_x, N] = size(x);
    end
    
    y = sum( x.* (A*x), 1);        

end