% function testTtest

    randn('state', 10)
    % random sample, with mean .1
    y = normrnd(.5, 1, 100,1);
    mu0 = 0;
    nbin = 10;
    figure(1)
    hist(y, nbin);
    m = mean(y);
    [h1,p1,ci1] = ttest(y,mu0);
    drawVerticalLine(m, 'color', 'k', 'linestyle', ':', 'linewidth', 2);
    drawVerticalLine(mu0, 'color', 'r', 'linestyle', '-');
    drawVerticalLine(ci1, 'color', 'g', 'linestyle', ':');
    title(sprintf('Two-tailed: h = %d, p = %.3f, ci = [%.3f, %.3f]', h1, p1, ci1(1), ci1(2)));
    xlim([-1 1]);
    
    figure(2)
    hist(y, nbin);
    [h2,p2,ci2] = ttest(y,mu0, [], 'left');
    drawVerticalLine(m, 'color', 'k', 'linestyle', ':', 'linewidth', 2);
    drawVerticalLine(mu0, 'color', 'r', 'linestyle', '-');
    drawVerticalLine(ci1, 'color', 'g', 'linestyle', ':');
    drawVerticalLine(ci2, 'color', 'b', 'linestyle', ':');
    title(sprintf('Left-tailed: h = %d, p = %.3f, ci = [%.3f, %.3f]', h2, p2, ci2(1), ci2(2)));
    xlim([-1 1]);
    
    figure(3)
    hist(y, nbin);
    [h3,p3,ci3] = ttest(y,mu0, [], 'right');
    drawVerticalLine(m, 'color', 'k', 'linestyle', ':', 'linewidth', 2);
    drawVerticalLine(mu0, 'color', 'r', 'linestyle', '-');
    drawVerticalLine(ci1, 'color', 'g', 'linestyle', ':');
    drawVerticalLine(ci3, 'color', 'm', 'linestyle', ':');
    title(sprintf('Right-tailed: h = %d, p = %.3f, ci = [%.3f, %.3f]', h3, p3, ci3(1), ci3(2)));
    xlim([-1 1]);
    
    
    
    
    
% end