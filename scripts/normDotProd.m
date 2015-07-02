function r = normDotProd(x,y)
    normx = norm(x(:));
    normy = norm(y(:));
    if (normx == 0) || (normy == 0)
        r = nan;        
    else
        x = x(:)/normx;
        y = y(:)/normy;
        r = dot(x,y);   
    end
    if isnan(r)
       3; 
    end
end