function ks_stat = lillieKS_stdCut(x, nStd)
    if nargin < 2
        nStd = 3;
    end
    warning('off', 'stats:lillietest:OutOfRangeP')
    x = x(:);
    m = mean(x);
    s = std(x);
    idx_include = (x > m - s*nStd) & (x < m + s*nStd);
    [h,p,ks_stat] = lillietest(x(idx_include));

end