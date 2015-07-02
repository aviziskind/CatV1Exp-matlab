function compareWvfms(Gid, clust1, clust2)

    
    S_wvfm = load(getFileName('mwvfm', Gid));
    if clust1 > clust2
        [clust1, clust2] = deal(clust2, clust1);
    end
    
    clustIds = [S_wvfm.meanWaveforms.cellId];
    wvfm_idx1 = find(clustIds == clust1, 1);
    wvfm_idx2 = find(clustIds == clust2, 1);    
    
    if isempty(wvfm_idx1) || isempty(wvfm_idx2)
        error('invalid cellid');
    end
    wvfmSclFactor = 1/4;
    figure(55); clf; 
%     subplot(2,1,1);
    t_ms = S_wvfm.t_ms; t_ms_ext = 1:(length(t_ms)*4);
    plot(t_ms_ext, [S_wvfm.meanWaveforms([wvfm_idx1 wvfm_idx2]).wvfm_raw]*wvfmSclFactor)
    xlim(t_ms_ext([1, end]));
    title(sprintf('Gid = %d. clustIds = %d, %d', Gid, clust1, clust2));
    
    return;
    
    S_ccg = load(getFileName('clusterData2', Gid));
    ccgData = S_ccg.data.ccgData;

    ccg_idx1 = find(S_ccg.data.clustIds == clust1, 1);
    ccg_idx2 = find(S_ccg.data.clustIds == clust2, 1);
    
    [allRanges_ms, allNbins] = getGlobals('isi_allRanges_ms', 'isi_allNbins');
    ccg_range_ms = 16;
    ccg_nbins = 80;
    rng_idx = find(allRanges_ms == ccg_range_ms, 1);
    nbin_idx = find(allNbins == ccg_nbins, 1);
    
    ccgBinVals = masterBin2Bin(ccgData.ccgs_mb{ccg_idx2, ccg_idx1}, ccgData.mbIdxs{rng_idx, nbin_idx});
    ccgBinEdges = linspace(-ccg_range_ms, ccg_range_ms, ccg_nbins+1);
    ccgBinCent = binEdge2cent(ccgBinEdges);
            
    subplot(2,1,2);
    bar(ccgBinCent, ccgBinVals, 1);
    xlim(ccgBinEdges([1, end]));
    

end