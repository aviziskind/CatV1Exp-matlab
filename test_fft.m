% 
% L = 200;
% w1 = 2;
% w2 = 20;
% f = rand(1, L)-.5;%cos(w*[0:L-1]*2*pi/L) + 2*cos(w2*[0:L-1]*2*pi/L);
% 
% 
% F = fft(f);
% figure(2); plot(f, 'o-');
% figure(3); plot(abs(F), 'o-'); 
% [fr, pw] = powerSpectrum(f);
% figure(4); plot(fr, pw, 'ro-');
% 

randn('seed', 2)
T = 8;
t1 = 0:T-1;
w = 1;
% f = cos(w*[0:L-1]*2*pi/L) ;
f1 = randn(T,1);

F1 = fft(f1);
[y, t2, f2] = fourierInterp(f1, 10, 2, 'harmonics');
[y, t3, f3] = fourierInterp(f1, 10, 2, 'spline');
% figure(1); clf; plot(t1, f1', 'bo-', t2, f2, 'rs:', t3, f3, 'k*:')
figure(1); clf; plot(t1, f1', 'bo-', t2, f2, 'rs:')

% [sum(abs(f2(1:2:end)-f1)), sum(abs(f3(1:2:end)-f1))]