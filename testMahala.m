% testMahala

N = 100;
x = randn(1,N);
y = randn(1,N)+2*x;

X = [x;y];

C = cov(X');
[V,D] = eig(C);

A = sqrt(D)\V';
Y = A*X;

figure(42);  clf; plot(X(1,:),X(2,:), 'b.'); hold on;
                  plot(Y(1,:),Y(2,:), 'g.')
figure(43);  clf; plot(Y(1,:),Y(2,:), 'g.')

r = randperm(N);
r = r(1:20);
[X1, Y1] = linesFromAtoB(X(1,r), X(2,r), Y(1,r), Y(2,r));
figure(42); hold on;
plot(X1, Y1, 'ro-');

