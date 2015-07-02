

clear all;
D = 2;

M1 = rand(D,1);
M2s = ones(D,1)*[-1.5:.5:1.5];

XX1 = randn(100,D);
XX2 = bsxfun(@plus, randn(100,D), M2s(:,1)');
C1 = cov(XX1);
C2 = cov(XX2);

figure(1); clf;
plot(XX1(:,1), XX1(:,2), 'b.', XX2(:,1), XX2(:,2), 'g.');


N = size(M2s,2);
ov_quad = zeros(1,N);
ov_anl = zeros(1,N);

gauss1 = @(X) gaussianN(X, M1, C1);
progressBar('init-', length(M2s))
for mi = 1:length(M2s)
    progressBar;
    M2 = M2s(:,mi);    
    nstd = 3;
    bnd = [M1-nstd*det(C1), M1+nstd*det(C1), M2-nstd*det(C2), M2+nstd*det(C2)];
    bnd1 = min(bnd,[], 2);
    bnd2 = max(bnd,[], 2);
    
    
    gauss2 = @(X) gaussianN(X, M2, C2);
    switch D
        case 1,  q_calc = quad(@(x) gauss1(x).*gauss2(x), bnd1(1), bnd2(1));
        case 2,  
            f1 = @(x,y) vectorFuncHandle(gauss1, x,y);
            f2 = @(x,y) vectorFuncHandle(gauss2, x,y);            
            q_calc = dblquad(@(x,y) f1(x,y).*f2(x,y), bnd1(1), bnd2(1), bnd1(2), bnd2(2));
        case 3,  
            f1 = @(x,y,z) vectorFuncHandle(gauss1, x,y,z);
            f2 = @(x,y,z) vectorFuncHandle(gauss2, x,y,z);            
            q_calc = triplequad(@(x,y,z) f1(x,y,z).*f2(x,y,z), bnd1(1), bnd2(1), bnd1(2), bnd2(2), bnd1(3), bnd2(3));
            
    end
    ov_quad(mi) = q_calc;
    
%     switch D
%         case 0,
%             s1 = sqrt(C1); s2 = sqrt(C2);
%             Z = 1/(2*pi*s1*s2);
%             A = 1/(2*s1^2) + 1/(2*s2^2);
%             B = M1/(s1^2) + M2/(s2^2);
%             C = -M1^2/(2*s1^2) - M2^2/(2*s2^2);
%             
%             q_an = Z*sqrt(pi/A)*exp(B^2/(4*A)+C);
%             ov_anl(mi) = q_an;
%         case {1,2,3}
%             C1i = inv(C1);
%             C2i = inv(C2);
% 
%             Z = 1/sqrt(((2*pi)^(2*D))*det(C1)*det(C2));
%             A = (C1i + C2i)/2;
%             B = (M1'*C1i + M2'*C2i)';
%             C = -(M1'*C1i*M1 + M2'*C2i*M2)/2;
% 
%             q_an =  Z * exp(.25*B'*inv(A)'*B+C) * sqrt(pi^length(A)*det(inv(A)));
%     end
    ov_anl(mi) = quadProdGaussians(M1, C1, M2, C2);     
            
end

figure(10); clf;
plot(M2s, ov_quad, 'bo-', M2s, ov_anl, 'r.:');

figure(11); clf;
plot(ov_quad, ov_anl, 'r.');
