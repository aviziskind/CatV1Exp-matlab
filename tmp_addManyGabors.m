%%
N = 10000;  
x = -3:.01:3;
nx = length(x);
ys = zeros(N, nx);
% phis = linspace(0, 360, N);
phis = rand(N, 1)*360;

for i = 1:N
    ys(i,:) = gabor1D(1, 0, 2, 3, phis(i) , x);    
end
% figure(1); clf;
% plot(ys');

%%

C = cumsum(ys, 1);

figure(2); clf; hold on;
plot(1:N, mean(C, 2), 'b-');
figure(3); clf; hold on;
plot(1:N, max(C,[], 2), 'r-');

%%
figure(4); clf; hold on;
nCol = 256;
j = jet(nCol);
for i = 1:2:N
    plot(sum(ys(1:i,:, 1), 1)', 'color', j(mod(i-1, nCol)+1,:))
    drawnow;
end



% ns = round(logspace(1, 3, N)); y = arrayfun(nrm, ns); plot(ns, y, '.-')