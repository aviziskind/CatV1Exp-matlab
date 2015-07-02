n40 = 178;
n60 = 1274;

n4 = 199;
n8 = 272;

Nbins = 7;

% [dphis, P] = deltaPhiNull(nPh, 6000);
% [dphis, P] = deltaPhiNull([4 8], [n4, n8]);
[dphis, P] = deltaPhiNull([40 60], [n40, n60]);
P = round(P*100);


if ~all(P == round(P))
    fprintf('not round numbers\n');
    return;
end

L = (180/(Nbins-1))/2;
binE= linspace(0-L, 180+L, Nbins+1);

dp = arrayfun(@(dph, p) repmat(dph, [1, p]), dphis, P, 'un', 0);
dp = [dp{:}];
figure(3);
binC = binEdge2cent(binE);
n = histcnt(dp, binE);
bar(binC,n);
set(gca, 'xtick', round(binC))
xlim([binE(1), binE(end)]);
title(num2str(mean(dp)))

% 40:  5!      9?  11?
% 60:  5?  7!  9?  11!  13?  17?