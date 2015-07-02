% function testCalcs

n = 1000000;
theta1 = rand(1, n)*pi; 
theta2 = rand(1, n)*pi;
phi1 = rand(1,n)*2*pi;
phi2 = rand(1,n)*2*pi;
r1 = 


v1 = randn(3,n);
[theta, phi, rho] = cart2sph(v1(1,:), v1(2,:), v1(3,:))


v2 = randn(3,n);
ccs = pearsonR_v(v1, v2);



d2 = cos(phi1-phi2);

d = cos(phi1).*sin(theta1).*cos(phi2).*sin(theta2) + sin(phi1).*sin(theta1).*sin(phi2).*sin(theta2) + cos(theta1).*cos(theta2);




v1 = randn(3,10000);
v1m = bsxfun(@minus, v1, mean(v1, 1));

a = [1, 1, 1]';
b = [1, -1, -1]';
c = [1, -1,  1]';

P = GramSchmidt([a,b,c]);
v1p = P'*v1m;

v1r = v1p(2:3,:);
ths = atan2(v1r(1,:), v1r(2,:));
rs = normV(v1r,1).^2;

figure(1); hist(ths, 30)
figure(2); normhist(rs, 30);

rs = normV(v1,1).^2;
v1

3;
return;





n = 200;
as = [linspace(0, .8, n/2), linspace(.8, 1, n/2)];


mus = zeros(1, n);
sigs = zeros(1,n);


th1 = linspace(0, 2*pi, 100); th1= th1(1:end-1);
th2 = linspace(0, 2*pi, 100); th2= th2(1:end-1);

[all_th1, all_th2] = meshgrid(th1, th2);
all_th1 = all_th1(:);
all_th2 = all_th2(:);

v1 = [cos(all_th1)'; sin(all_th1)'];
v2 = [cos(all_th2)'; sin(all_th2)'];

vb = [1;0];

progressBar('init-', length(as))
for i = 1:length(as)
    progressBar;

    a = as(i);
    
    v1p = bsxfun(@plus, a*vb, (1-a)*v1);
    v2p = bsxfun(@plus, a*vb, (1-a)*v2);

    
%     v1p = bsxfun(@rdivide, v1p, normV(v1p, 1));
%     v2p = bsxfun(@rdivide, v2p, normV(v2p, 1));
% 
%     ccs = sum(v1p .* v2p, 1);

    ccs = pearsonR_v(v1p, v2p);
    
    
    binE = linspace(-1, 1, 30); binC = binEdge2cent(binE);
    binV = histcnt(ccs, binE);
    figure(1); clf; bar(binC, binV, 1);

%     fprintf('mean = %.2f   var = %.2f\n', mean(ccs), var(ccs))

    mus(i) = mean(ccs);
    sigs(i) = var(ccs);
    
    
end


figure(2); clf;
plot(mus, sigs, 'g.');
% vb = 1/sqrt(2)*[1;1];

% a = .3;
% ccs_b1 = a*vb + (1-a)*


% v1 = [xg, yg];
% 
% v1 = 



79-
99-



