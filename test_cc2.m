function test_cc2

    n = 8;
    B = 1000;
    f1odcs = [.1:1.9];
    nF1 = length(f1odcs);
    eta = .1;



    ccs = zeros(1,nF1);
    x_end = 1.6;
    dy = [.2:.01:x_end];
    f1s = zeros(size(dy));

    for fi = 1:length(dy)
        ph = linspace(0, 2*pi, n+1); ph = ph(1:n);
        r1 = gaussian(1:n, n/2, dy(fi));
    %     r1 = rectified(sin(ph)+dy(fi) );
        f1s(fi) = getF1oDC(ph, r1, 2*pi);    
    end


    figure(545); clf;
    plot(dy, f1s, '.'); hold on;
    % hyp = @(beta, x) beta./abs(x+1e-5);
    % beta = nlinfit(dy, f1s, hyp, 1);
    % fplot(@(x) hyp(1, x), [.5 10], 'r');
    % xlim([0 20])
    % gs = @(beta, x) beta(1) * gaussian(x, 0, beta(2));
    p = polyfit(dy, f1s, 1);
    p_inv = [1/p(1), -p(2)/p(1)]
    % beta = nlinfit(dy, f1s, gs, [40, 8]);
    fplot(@(x) polyval(p, x), [.2 x_end], 'r');
    xlim([0 x_end])
    %}





        n = 8;
        ph = linspace(0, 2*pi, n+1); ph = ph(1:n);

        f1odcs1 = [0:.05:2];
        f1odcs2 = zeros(size(f1odcs1));
        for i = 1:length(f1odcs1);
            tc = getTCwithF1oDC(f1odcs1(i), 8);
            f1odcs2(i) = getF1oDC(ph, tc, 2*pi);

        end
        figure(434); plot(f1odcs1, f1odcs2, 'o:')
end



function tc = getTCwithF1oDC(f1odc, n)
    ph = linspace(0, 2*pi, n+1); ph = ph(1:n);

    if f1odc <= 1  % use sinusoid + const (hyperbolic)
        cnst = 1/(1e-5+abs(f1odc));
        tc = rectified(sin(ph)+ cnst);
    elseif f1odc > 1 % use gaussian with changing sigma (linear)
        p_inv = [-1.2279    2.7843];
        sig = polyval(p_inv, f1odc);
        tc = gaussian(1:n, n/2, sig);
    end
end
