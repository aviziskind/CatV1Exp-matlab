function testVarianceOfCosGaussian



    M = 10;
    S = 3;
    lim = 15;
    dx = .01;
    figure(55); clf;
    f = @(x) gaussian(x, M, S);
    fcos = @(x) (f(x) .* cos(10*x)).^2;

    xs = M-lim:dx:M+lim;
    ys1 = f(xs);
    ys2 = fcos(xs);

    mn1 = distribMean(xs, ys1) %#ok<NOPRT>
    mn2 = distribMean(xs, ys2) %#ok<NOPRT>

    figure(4); clf(4);
    plot(  (ys1*dx) .* xs , 'b');
    hold on;
    plot(  (ys2*dx) .* xs , 'g');
    

	std1 = distribStd(dx, ys1)
    std2 = distribStd(dx, ys2)
    figure(2); clf(2);
    plot(xs, ys1/sum(ys1*dx), 'o-');
    hold on;
    plot(xs, ys2/sum(ys2*dx), 'g.-');



end