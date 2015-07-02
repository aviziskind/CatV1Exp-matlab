

M1 = [0;0];
M2s = [1;1]*[-1.5:.5:1.5];

C1 = cov(randn(100,2));
C2 = cov(randn(100,2));

% 
% C1 = [5, .1; .3, 4];
% C2 = [2, .3; .1, 3];
% 

N = size(M2s,2);
% ov_min = zeros(1,N);
ov_quad = zeros(1,N);
ov_anl = zeros(1,N);

gauss1 = @(X) gaussianN(X, M1, C1);


f1 = @(x,y) vectorFuncHandle(gauss1, x,y);

fmesh(f1, [-3, 3], [-3, 3])

k = 2;

for mi = 1:length(M2s)
    M2 = M2s(:,mi);
%     ov_min(mi) = gaussiansOverlap(m1, s1, m2, s2);
    gauss2 = @(X) gaussianN(X, M2, C2);

    f1 = @(x,y) vectorFuncHandle(gauss1, x,y);
    f2 = @(x,y) vectorFuncHandle(gauss2, x,y);
    f_prod = @(x,y) f1(x,y).*f2(x,y);
    q_calc = dblquad(f_prod, -10, 15, -10, 15);

    
%     tst1 = dblquad( exp( 
    
    ov_quad(mi) = q_calc;

    C1i = inv(C1);
    C2i = inv(C2);

    Z = 1/sqrt(((2*pi)^(2*k))*det(C1)*det(C2));
    A = (C1i + C2i)/2;
    B = (M1'*C1i + M2'*C2i)';
    C = -(M1'*C1i*M1 + M2'*C2i*M2)/2;
    
    %     
%     A = 1/(2*pi*s1*s2);
%     B = 1/(2*s1^2) + 1/(2*s2^2);
%     C = m1/(s1^2) + m2/(s2^2);
%     F = -m1^2/(2*s1^2) - m2^2/(2*s2^2);
    
    q_an =  Z * exp(.25*B'*inv(A)'*B+C) * sqrt(pi^length(A)*det(inv(A)));

%     q_an = Z*sqrt( ((2*pi)^k) * det(inv(A)) )*exp(C + B'*inv(A)'*B );
    
    
    ov_anl(mi) = q_an;
end

figure(1); clf;
plot(M2s, ov_quad, 'bo-', M2s, ov_anl, 'r.-');
set(gca, 'yscale', 'log')

% figure(2)
% plot(ov_quad, ov_anl, 'r.');
