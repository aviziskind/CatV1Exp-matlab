function s = str_moments(x, n)

    fmt = '%.2f';
    
    if nargin < 2
        n = 2;
    end
    
    m = zeros(n,1);
    m(1) = nanmean(x);
    N = nnz(~isnan(x));
    stdx = nanstd(x); 
    sem = stdx/sqrt(N);
    if n >= 2
        m(2) = stdx;
    end
    if n >= 3
        m(3) = skewness(x);
    end
    if n >= 4
        m(4) = kurtosis(x);
    end
        
    s = ['[' num2str(m(1), fmt) ' (\pm ' num2str(sem, fmt) '), '];
    for i = 2:n-1
        s = [s, sprintf([fmt ', '], m(i))]; %#ok<AGROW>
    end
    s = [s sprintf([fmt ']'], m(n))];
        
end