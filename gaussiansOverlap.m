function s = gaussiansOverlap(m1, s1, m2, s2, showFlag)

    showWorking = exist('showFlag', 'var') && ~isempty(showFlag);
    
    methodId = 1;

    if m1 > m2
        [m1, s1, m2, s2] = deal(m2, s2, m1, s1); % swap.
    end
    
    nStd = 5;
    xs_extrm = [m1-s1*nStd, m1+s1*nStd, m2-s2*nStd, m2+s2*nStd];
    
    x_bnd = [min(xs_extrm), max(xs_extrm)];
    
    xs = linspace(x_bnd(1), x_bnd(2), 100);
    y1 = gaussian(xs, m1, s1);
    y2 = gaussian(xs, m2, s2);

    f1 = @(x) gaussian(x, m1, s1);
    f2 = @(x) gaussian(x, m2, s2);
    
    if (methodId == 1)
        f_min = @(x) min(f1(x), f2(x));
        s = quad(f_min, x_bnd(1), x_bnd(2));
    end
    if (methodId == 2) || showWorking
        diff_fs = @(x) f1(x)-f2(x);
        x_int = fzero(diff_fs, (m1+m2)/2);
        assert(ibetween(x_int, m1, m2));
        erfc1 = @(x) .5*erfc(x/sqrt(2)); % integral of gaussian with variance 1 (instead of 1/2)        
        area_half1 = erfc1(abs(x_int-m2)/s2);
        area_half2 = erfc1(abs(x_int-m1)/s1);
        s = area_half1 + area_half2;    
    end
    
    
    if showWorking
        figure(17); clf; hold on;
        plot(xs, y1, 'bo', xs, y2, 'go')
        fplot(f1, x_bnd, 'b');
        fplot(f2, x_bnd, 'g');
        fplot(f_min, x_bnd, 'r');
        plot(x_int, f1(x_int), 'rs');        
    end    

end