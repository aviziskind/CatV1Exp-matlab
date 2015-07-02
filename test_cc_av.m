n = size(ccs_smp,3);
meanDiffs = zeros(1, n);
cc_av = mean(ccs_smp,3);
for i = 1:n
    Ci = abs( mean(ccs_smp(:,:,1:i), 3) - cc_av );
    meanDiffs(i) = sum(Ci(:));    
end
figure(505); clf;
plot(1:n, meanDiffs);
