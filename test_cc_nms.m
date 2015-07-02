function test_cc_nms
% test relationship between cc and non-mean subtracted cc


mu = 1;

n = 8;
B = 2000;

[ccs, ccs_nms, ccs_dist] = deal(zeros(1,B));
ph = linspace(0, 360, n+1); ph=ph(1:n);


doTests = [2];

if any(doTests == 1)
    tic;
    for b = 1:B

    %     tc1a = rand(n, 1);
    %     tc1b = rand(n, 1);
        tc1a = exprnd(mu, n, 1);
        tc1b = exprnd(mu, n, 1);
        tc1a = tc1a/norm(tc1a);
        tc1b = tc1b/norm(tc1b);
        cc = corr(tc1a(:), tc1b(:));    
        cc_nms = corr_nms(tc1a(:), tc1b(:));
    %     if cc < -.75 && cc_nms > .5
        if cc > 0 && cc_nms > .5
            figure(567); clf;
            plot(ph, tc1a, 'bo-', ph, tc1b, 'go-');
            set(gca, 'xtick', ph(1:2:end))
            title(sprintf('cc = %.2f. dot = %.2f', cc, cc_nms))
            3;
        end
        ccs(b) = cc;
        ccs_nms(b) = cc_nms; 
        ccs_dist(b) = corr_dist(tc1a(:), tc1b(:));

    end


    figure(424); plot(ccs, ccs_nms, 'b.'); xlabel('cc'); ylabel('dot');
    figure(425); plot(ccs, ccs_dist, 'g.'); xlabel('cc'); ylabel('dist');
    figure(426); plot(ccs_nms, ccs_dist, 'r.'); xlabel('dot'); ylabel('dist');
    nbin = 20;
    figure(445); hist(ccs, nbin); title(sprintf('cc (n = %d): %.2f \\pm %.2f', n, mean(ccs), stderr(ccs)));
    figure(446); hist(ccs_nms, nbin); title(sprintf('dot (n = %d): %.2f \\pm %.2f', n, mean(ccs_nms), stderr(ccs_nms)));
    figure(447); hist(ccs_dist, nbin); title(sprintf('dist (n = %d): %.2f \\pm %.2f', n, mean(ccs_dist), stderr(ccs_dist)));
end

if any(doTests == 2)  % show how adding noise affects tuning curves of various f1/dcs
    
    etas = [.05:.15:.7];
    f1odcs = [0:.1:2];
    nE = length(etas);
    nF = length(f1odcs);
    nT = 50;
    n = 8;
    
    ccs = zeros(nE, nF, nT);
    dots = zeros(nE, nF, nT);
    ph = linspace(0, 360, n+1); ph=ph(1:n);

    for ei = 1:length(etas)
        for fi = 1:length(f1odcs)
            for ti = 1:nT
                tc = getTCwithF1oDC(f1odcs(fi), n);            
                tc1 = rectified( tc + randn(1,n)*mean(tc)*etas(ei) );
                tc2 = rectified( tc + randn(1,n)*mean(tc)*etas(ei) );
                
                cc = doPearsonCorr(tc1, tc2);
                dt = normDotProd(tc1, tc2);
                ccs(ei, fi, ti) = cc;
                dots(ei, fi, ti) = dt;        
                
                if (ti == 1) && (ei == 3) && (fi == 3)
                    figure(888);                    
                    plot(ph, tc, 'k.:', ph, tc1, 'bo-', ph, tc2, 'go-');
                    set(gca, 'xtick', ph(1:2:end))
                    s1 = sprintf('F1/DC = %.2f. \\eta = %.2f', f1odcs(fi), etas(ei));
                    s2 = sprintf('cc = %.2f. dot = %.2f', cc, dt);
                    ylim([0 inf]);
                    title({s1, s2});
                    3;
                end
                
            end
        end
    end
    
    figure(879); clf; hold on
    for ei = 1:length(etas)
        y = mean(ccs(ei, :, :), 3);
        e = std(ccs(ei, :, :), [], 3);
        errorbar(f1odcs, y, e, ['-o' color_s(ei)]);        
    end
    box on;
    xlabel('F1/DC'); ylabel('CC');
    ylim([-1 1]);
    xlim([-.1, 2.1]);
    l_array = legendarray('\eta = ', etas);
    legend(l_array, 'location', 'SE');
    
    figure(880); clf; hold on
    for ei = 1:length(etas)
        y = mean(dots(ei, :, :), 3);
        e = std(dots(ei, :, :), [], 3);
        errorbar(f1odcs, y, e, ['-s' color_s(ei)]);        
    end
    box on;
    ylim([0 1]);
    xlim([-.1, 2.1]);
    legend(l_array, 'location', 'SE');
    xlabel('F1/DC'); ylabel('Dot products');
    
end
    
    

end

function r = corr_nms(x, y)
    x = x/norm(x);
    y = y/norm(y);
    r = dot(x,y);                
end

function r = corr_dist(x, y)
    x = x/norm(x);
    y = y/norm(y);
    r = norm(x-y);
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
    tc = tc/norm(tc);
end

