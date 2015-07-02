Tmax = 100;

Twind = 10;
dt = .1;
tau = 4;
eta = 10;

t_full = 0:dt:Tmax;
t_wind = 0:dt:Twind;
nmax = length(t_wind);

y = zeros(1, length(t_full));
y([1 2]) = [1 2];
for i = 3:length(t_full)
    n = min(i-1, tau);
    p = polyfit(1:n, y(i-n:i-1)-y(i-n), 1);    
    y(i) = polyval(p, n+1)+randn*eta;
     
end

c = xcorr_wrapped(y, y, nmax);
figure(1); clf;
plot(t_full, y);
figure(2); clf;
plot(t_wind, c, '-o');