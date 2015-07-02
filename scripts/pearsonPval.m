function p = pearsonPval(tail, rho, n)
    %PVALPEARSON Tail probability for Pearson's linear correlation.
    % copied from "corr.m"
    t = sign(rho) .* Inf;
    k = (abs(rho) < 1);
    t(k) = rho(k).*sqrt((n-2)./(1-rho(k).^2));
    switch tail
    case 'b' % 'both or 'ne'
        p = 2*tcdf(-abs(t),n-2);
    case 'r' % 'right' or 'gt'
        p = tcdf(-t,n-2);
    case 'l' % 'left' or 'lt'
        p = tcdf(t,n-2);
    end
    
end