% testSaddlePoint

N = 20;
xs = 1:N;
ys = 1:N;
[xs_grid, ys_grid] = meshgrid(xs, ys);

G1 = reshape(gaussianN([xs_grid(:)'; ys_grid(:)'], [5;10], 5*eye(2)), N, N);
G2 = reshape(gaussianN([xs_grid(:)'; ys_grid(:)'], [14;10], 4*eye(2)), N, N);

G = G1 + G2;

dF1 = diff(G, 2, 1);
dF2 = diff(G, 2, 2);
dF1 = dF1(:,2:N-1);
dF2 = dF2(2:N-1,:);
% dF_mag = normV(cat(3, dF1, dF2), 3);
dF_mag = dF1 .* dF2;
% diff(G, 1, 2);
figure(1); surf(G);
figure(2); surf(dF1);
figure(3); surf(dF2);
figure(4); surf(dF_mag);
