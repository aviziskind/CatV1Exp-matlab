function [meanMID_cc, allMID_ccs] = getMeanMIDcc(v_MID)
    [m,n,nJacks] = size(v_MID);
    MID_vecs = reshape(v_MID, [m*n, nJacks]);
    mid_CCs = pearsonRm(MID_vecs);
    idx_lower = find( tril(ones(nJacks), -1) );
    
    allMID_ccs = mid_CCs(idx_lower); %#ok<FNDSB>
    meanMID_cc = mean(allMID_ccs);

end