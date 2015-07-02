function D = invDiagThreshold(D, th)
    n = length(D);
    diagInds = sub2ind([n n], 1:n, 1:n);
    
    diagIndsAboveTh = diagInds ( abs(diag(D)) > th );   
    D(diagIndsAboveTh) = 1./D(diagIndsAboveTh);
    
%     diagIndsBelowTh = diagInds ( abs(diag(D)) <= th );   
%     D(diagIndsBelowTh) = 1;
end