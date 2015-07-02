function [n, eff] = matchNumbersToRatios(n, r)
% x_ratios: [.7  .3 ]
% y_vals:   [90  10 ]
L = length(n);
n_start = n;

pairs = nchoosek(1:L,2);
pairs = [pairs; fliplr(pairs)];

for i = 1:size(pairs,1);
    i1 = pairs(i,1); i2 = pairs(i,2);
    
    if (n(i1)-1)/n(i2) > r(i1)/r(i2)   % n(1) is too high
       n(i1) = (r(i1)/r(i2))*n(i2);
    end   
end
% assert(all(r/sum(r) == n/sum(n)))
assert( sum(abs(r/sum(r) - n/sum(n))) < 1e-10)
n = round(n);

eff = sum(n)/sum(n_start);

end