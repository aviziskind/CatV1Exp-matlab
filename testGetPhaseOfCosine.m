

xs = -5:.05:6;
k = .1 + rand*3;
x0 = randn;
phi = rand*2*pi;
ys = cos(k*(xs-x0) + phi); 
phi2 = getPhaseOfCosine(xs, ys, k, x0);
ys2 = cos(k*(xs-x0) + phi2); 

figure(1); clf
plot(xs, ys, 'b.-', xs, ys2, 'g.-');
title(sprintf('k = %3.2g, x0 = %3.2g, phi_1 = %3.3g, phi_2 = %3.3g', k, x0, phi, phi2));

