function testSparse

L = 120;
N = 3690;

m = randi(L, N, 1);
n = randi(L, N, 1);
v = randn(N, 1);

A_s = sparse(m, n, v, L, L);
A_f = full(A_s);


% test retrieval time
tic;
for i = 1:N
    A_s(m(i), n(i));
end
toc;
t1 = toc;

tic;
for i = 1:N
    A_f(m(i), n(i));
end
toc;
t2 = toc;
% t1/t2

% test find time
tic;
for i = 1:N;
    [r,c,v] = find(A_s);
end
toc;
t1 = toc;


tic;
for i = 1:N
    [r,c,v] = find(A_f);
end
toc;
t1 = toc;
t1/t2