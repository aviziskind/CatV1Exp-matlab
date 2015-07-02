function h = stairs2(x,y, varargin)

    ext_extrp = @(x) [x(:); x(end)+diff(x([1:2]))];        
    ext_rep   = @(x) [x(:); x(end)];        
    
    h = stairs(ext_extrp(x), ext_rep(y), varargin{:});    
end
