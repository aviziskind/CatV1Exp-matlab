function Y = cat2cos(B, N1, N2)
%     [A0, A1, A2, A3, phi1, phi2, phi3, B0, B1, B2, B3, psi1, psi2, psi3]
    
    N = length(B);
%     N = 14;
    n_each = N/2;
    m = (n_each-1)/2;

%     Nb = length(B);
%     m = Nb/2 - 1;
%     n = m/2+1;
    B1 = B(1:n_each);
    B2 = B(n_each+1:N);
   
    Amps1 = B1(1:m+1);
    Amps2 = B2(1:m+1);
    Phis1 = [0, B1(m+2:n_each)];
    Phis2 = [0, B2(m+2:n_each)];
    
%     A1   = B(    1 :   m);
%     phi1 = B(  m+1 : 2*m);
%     A2   = B(2*m+1 : 3*m);
%     phi2 = B(3*m+1 :   N);
    
    X1 = (0:N1-1)*2*pi/N1;
    X2 = (0:N2-1)*2*pi/N2;
    
    Y = [ones(1, N1)*Amps1(1), ones(1, N2)*Amps2(1)];
    for i = 1:m
        Y = Y + [Amps1(1+i)* cos(i*X1+Phis1(i+1)), Amps2(1+i)* cos(i*X2+Phis2(i+1))];
    
    end
    Y = Y/max(abs(Y));
end
% A*cos(idx1*2*pi/Nt1 + phi);