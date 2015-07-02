




T = 8;
t1 = 0:T-1;
% t2 = 0:.5:T-.5;
t1ext = 0:2*T-1;

% f = cos(w*[0:L-1]*2*pi/L) ;
f1 = randn(T,1);

[y, t2, f2] = fourierInterp(f1, 10, 2, 'spline');

F1 = fft(f1); itpM = 2;
mid = round(T/2);
F_itp = [F1(1:mid); zeros(length(F1)*(itpM-1), 1); F1(mid+1:end)]*itpM/2;


figure(1); clf; plot(t1, f1', 'bo-', t2, f2, 'rs:')
figure(2); clf; plot(t1ext, fft(f2)/2, 'rs:', t1ext, F_itp, 'b.-')
drawVerticalLine([mid, mid-1+length(F1)*(itpM-1)], 'linestyle', ':')
zeroaxes;
% [y, t3, f3] = fourierItp(f1, 10, 'spline', 2);
% figure(1); clf; plot(t1, f1', 'bo-', t2, f2, 'rs:', t3, f3, 'k*:')
