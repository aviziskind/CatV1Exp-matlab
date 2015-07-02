function n_wgt = weightedBinCounts(bin, wgts, N)

    idx_ok = find( ibetween(bin(:)', 1, N) );

    n_wgt = zeros(1,N);
    for i = idx_ok
        n_wgt(bin(i)) = n_wgt(bin(i)) + wgts(i);
    end
end