function testRankSum

    
    x = randn(1000,1);
    y0 = randn(2000,1);
    
    L = min([x(:); y0(:)]) - 2;
    R = max([x(:); y0(:)]) + 2;
    
    bin_x = linspace(L,R,30);
    binC = binEdge2cent(bin_x);
    nx = histcnt(x, bin_x);
    ny = histcnt(y0, bin_x);
    
    
    figure(1);
    subplot(2,1,1);
    hx = bar(binC, nx);
    ht = title(' ');
    subplot(2,1,2);
    hy = bar(binC, ny);
    
    function showPlots(dy)
        y = y0 + dy;
        ny = histcnt(y, bin_x);
        set(hy, 'ydata', ny);
        
        [h, pt_0] = ttest2(x,y, [], 'both');
        [h, pt_L] = ttest2(x,y, [], 'left');
        [h, pt_R] = ttest2(x,y, [], 'right');    
        str_t = sprintf('t-test: 2-tailed: %.2g. Left: %.2g, Right: %.2g', pt_0, pt_L, pt_R);

        pW_0 = myRanksum(x,y, 'both');
        pW_L = myRanksum(x,y, 'left');
        pW_R = myRanksum(x,y, 'right');    
        str_U = sprintf('U-test: 2-tailed: %.2g. Left: %.2g, Right: %.2g', pW_0, pW_L, pW_R);        

        set(ht, 'string', {str_t, str_U})
    end


    manipulate(@showPlots, {'dy', [-2, 2], 0}, 'FigId', 2);


end