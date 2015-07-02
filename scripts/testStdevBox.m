function testStdevBox
    figure(1);
    clf;
    axis([-1 3, 0 5]);
    x1 = randn(100, 1);
    x2 = randn(100, 1) + 1;
    
    
    h1 = stdevBox(0, 1, 1, x1, .3, .4, 'r');
    h2 = stdevBox(0, 1, 2, x2, .3, .4, 'r');
%     stdevBox(x,y, stdev, sem, w1, w2, col, hv)

    pause(.5)
    x1 = x1+.2;
    x2 = x2-.1;
    h1 = stdevBox(h1, 1, 1, x1, .3, .4, 'r');
    h2 = stdevBox(h2, 1, 2, x2, .3, .4, 'r');

    
end