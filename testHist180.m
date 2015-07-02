x = (randi(5,1,100)-1)*45;
L = 45/4;
binEdges = linspace(-L, 180+L, 10);
[n, whichBins] = histcnt(x, binEdges);
binCenters = binEdge2cent(binEdges);
bar(binCenters, n, 1);
% set(gca, 'xtick', 0:45:180)
xlim([-L, 180+L])
