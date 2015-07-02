N = 200;
x = linspace(0,1,N)';
y = x * 0.2 + rand(N,1);
figure(1); clf;
plot(x,y, '.'); hold on;

p = polyfit(x,y, 1);
disp(['slope : ' num2str(p(1))]);
f = polyval(p, x);
plot(x, f, 'g-');

X = [ones(size(x)), x];
alpha = .01;
[b, bint] = regress(y, X, alpha)
