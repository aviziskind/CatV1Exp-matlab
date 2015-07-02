function b = pca_slope(x,y)
    x0 = mean(x);
    y0 = mean(y);
    xt = x-x0;
    yt = y-y0;
    
    v = [xt; yt];
    
    xx = dot(xt,xt);
    yy = dot(yt,yt);
    xy = dot(xt,yt);
    A = [xx, xy; xy, yy];
    
    [V,d] = eig(A);
    max_eig_id = indmax(diag(d));
   
    v = V(:,max_eig_id);
    b1 = v(2)/v(1);
    b0 = y0 - b1*x0;
    
    b = [b0, b1];
end