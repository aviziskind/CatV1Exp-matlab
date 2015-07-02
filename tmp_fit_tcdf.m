[8.8000   18.9300]
p = [-0.00385, 0.2135, 0.248, 0.303];
pval_tmp = @(tstat) 10.^(-polyval(p, tstat));

figure(545); clf; hold on;

t = [0:.2:22];
p = tcdf(-t,360-2);
plot(t,-log10(p), 'b.')

idx = (p==0);
p(idx) = arrayfun(pval_tmp, t(idx));
plot(t,-log10(p), 'ro')

return;


% x = x(~isinf(y));
% y = y(~isinf(y));



% 
% 
% x = [0:.01:24];
% y = arrayfun(@(t)  -log10( tcdf(-t,360-2) ) , x);
% x = x(~isinf(y));
% y = y(~isinf(y));



x2 = [0:.5:24];
% poly fit
% p = polyfit(x, y, 3)
p = [-0.00385, 0.2135, 0.248, 0.303];
f = polyval(p, x);
figure(545); clf; hold on;
plot(x,y, 'b.')
plot(x2,f, 'r:')


p = [-0.00385, 0.2135, 0.248, 0.303];
pval = @(tstat) 10.^(-polyval(p, tstat));





% logistic fit
% func = @(b, x)   b(4) * ( abs(x.^b(1)) ./ abs(x.^(b(2) + b(3))) );
% B = nlinfit(x, y, func, rand(4,1));
% figure(545); clf; hold on;
% plot(x,y, '.')
% plot(x, func(B, x), 'r:')
% 
% x2 = 1:1:30;


% B = nlinfit(x, y, func, [1, 1]);
% 
% func = @(b, x)   -log10( normcdf(-x, b(1), b(2) ) );
% B = [1 1];
% figure(545); clf; hold on;
% plot(x,y, '.')
% plot(x, func(B, x), 'r.')
% 
% 
% y = arrayfun(@(t)  -log10( tcdf(-t,360-2) ) , x);
% 
% P = normcdf(X,mu,sigma)