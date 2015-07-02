L = 5; N = 1e5;
% X = rand(L,N);
% Y = rand(L,N);
X = random('poiss', 1, L,N);
Y = random('poiss', 1, L,N);
ccs = pearsonR_v(X, Y);
binE = linspace(-1, 1, 51);
binC = binEdge2cent(binE);
binV = histcnt(ccs, binE);
bar(binC, binV, 1);

fprintf('Ndof = %d   1/var = %.3f\n', L-1, 1/var(ccs));