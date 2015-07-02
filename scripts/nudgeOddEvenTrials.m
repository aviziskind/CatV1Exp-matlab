function [r1, r2] = nudgeOddEvenTrials(r1,r2)
    r1 = r1(:); r2 = r2(:);
    [urows, nReps] = uniqueCount([r1, r2], 'rows');

    D = 1/50;
    if all(nReps == 1)
        return;
    end
    
    d = min(diff(r1));
    sz = size(r1);
    r1 = r1 + (d*D)*randn(sz);
    r2 = r2 + (d*D)*randn(sz);

end