function y = gaussian2D_hlp(x,y,M,C)
    x = x(:)';
    y = y(:)';
    nX = length(x);
    nY = length(y);
    if (nX == 1) && (nY == 1)
        y = gaussian2D([x;y], M,C);
    elseif (nX > 1) && (nY == 1)
        y = gaussian2D([x; ones(1,nX)*y], M,C);
    elseif (nX == 1) && (nY > 1)
        y = gaussian2D([x*ones(1,nY); y], M,C);
    elseif (nX == nY)
        y = gaussian2D([x;y], M,C);
    else
        error('mismatch of dimensions');
    end


end