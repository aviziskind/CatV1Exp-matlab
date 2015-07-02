X = randn(200,1)*5+10;
nStd = 1.5;
m = mean(X);
s = std(X);

figure(1);
clf; normhist(X, 30);
hold on
fplot(@(x) gaussian(x, m, s), xlim, 'r:')

drawVerticalLine([m + s*nStd, m - s*nStd], 'linestyle', '--', 'color', 'g')

X2 = abs(X-m)/s;
erfc1 = @(x) erfc(x/sqrt(2)); 

p = arrayfun(@(x) erfc1(x), X2);

p_th = erfc1(nStd);

idx1 = find( ~between(X, m-s*nStd, m+s*nStd) );
idx2 = find( p < p_th );
figure(2); clf;
plot( abs(X-m), p, 'b.')
isequal(idx1, idx2)
