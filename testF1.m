function testF1

    phis = 0;%:pi/10:2*pi;
    figure(1); clf; h1 = plot(0,0, 'b.'); hold on; h2 = plot(0,0, 'g-'); hold off;
    figure(3); clf; hold on; axis([-2 2, -2, 2]);
    for i = 1:length(phis)
        phi = phis(i);
        dt = .1;
        t = 0:dt:10*pi;
        T = 2*pi;

        w = 1;
        f = sin(w*t + phi) + randn(size(t));
        

        fourierTransform = @(w)   (1/sqrt(2))* sum( dt .* exp(1i * w * t) .* f) /T;
        fp = @(q) fourierTransform(q);
        fm = @(q) fourierTransform(-q);
        F = @(q) 2 * sqrt( (abs( fp(q) ))^2 + (abs( fm(q) ))^2 );    
%         F2 = @(q)  abs( fp(q) );
    
        figure(2); fplot(F, [-5 5]);
        figure(3);
        
        plot(real(fp(w)), imag(fp(w)), 'bo'); 
        plot(real(fm(w)), imag(fm(w)), 'ro'); 
        Q = -angle(fp(w)) + pi/2;
        
        set(h1, 'xdata', t, 'ydata', f);
        set(h2, 'xdata', t, 'ydata', sin(w*t + Q));
        
        drawnow;
%         pause(.01);
        
    end
        

    
end