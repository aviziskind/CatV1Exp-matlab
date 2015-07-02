% D = abs(randn(200,1));
% Dend = 1.3;
% D = D(D<=Dend);
% hist(D)

n1 = 20;
cdf1_flat = [linspace(0, Dend, n1)', linspace(0,1,n1)'];
[h,p] = kstest(D, cdf1_flat);
fprintf('p_ks = %.4g\n', p);

ns = [1e2, 1e4, 1e5, 1e6, 1e7];
for i = 1:length(ns)
    pdf1_flat = linspace(0,Dend, ns(i));
    [h,p] = ttest2(D, pdf1_flat);
    fprintf('p_t (n = %g) = %.4g\n', ns(i), p);
end

ns = [1e2, 1e4, 1e5, 1e6, 1e7];
for i = 1:length(ns)
    pdf1_flat = linspace(0,Dend, ns(i));
    p = ranksum(D, pdf1_flat);
    fprintf('p_U (n = %g) = %.4g\n', ns(i), p);
end

% n = 1e7;
% pdf1_flat = linspace(0,2, n);
% [h,p] = ttest2(D, pdf1_flat)