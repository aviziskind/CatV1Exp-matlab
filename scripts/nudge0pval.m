function p = nudge0pval(p, amt)
    if nargin < 2
        amt = 1e-50;
    end
    p(p == 0) = amt;
end