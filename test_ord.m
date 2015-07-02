for i = 1:100
shft = @(x,n) x([end-n+1:end, 1:end-n]);
N = 30;
x = randn(1,N);
y = randn(1,N);

ccs = arrayfun(@(n) pearsonR(x, shft(y,n)), 1:N);
dots = arrayfun(@(n) dot(x, shft(y,n)), 1:N);

idx1 = ord(ccs);
idx2 = ord(dots);
assert(isequal(idx1, idx2))
end