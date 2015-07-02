function [densHists, axesLims2D, axesTickLims2D] = getDensityHists(features, fetPairs, cellSpikeIdx, nBins, ignoreFirstClustFlag)
    nCells = length(cellSpikeIdx);

    nPairs = size(fetPairs,1);
    densHists.abs_lin(1:nPairs) = {zeros(nBins,nBins, nCells, 'single')};             
        
    % find x & y limits for each pair, for each co-ordinate basis
    if exist('ignoreFirstClustFlag', 'var') && ~isempty(ignoreFirstClustFlag) && (ignoreFirstClustFlag == true)
        clustIdxsForRange = 2:nCells;
    else
        clustIdxsForRange = 1:nCells;
    end
    dim_cat = find(size(cellSpikeIdx{1}) > 1, 1);
    allIdx = cat(dim_cat, cellSpikeIdx{clustIdxsForRange});
    fetLims = lims(features(allIdx,:), 0, 1);    
                     
    for pair_i = 1:nPairs
        [idx1, idx2] = deal(fetPairs(pair_i,1), fetPairs(pair_i,2));

        % find spike densities for each for each pair, for each co-ordinate basis                
        x_edges = linspace(fetLims(1, idx1), fetLims(2, idx1), nBins+1);
        y_edges = linspace(fetLims(1, idx2), fetLims(2, idx2), nBins+1);

        for cell_i = 1:nCells
            densityHist = hist2d([features(cellSpikeIdx{cell_i}, idx1), features(cellSpikeIdx{cell_i}, idx2)], x_edges, y_edges);
                        
            idx_nz = find(densityHist ~= 0);
            logDensityHist = zeros(size(densityHist));
            logDensityHist(idx_nz) = log10(densityHist(idx_nz));
                    
            densityHist_norm = densityHist / max( densityHist(:) );
            logDensityHist_norm = logDensityHist / max(logDensityHist(:) );

            densHists.abs_lin{pair_i}(:,:,cell_i) = densityHist;
            densHists.abs_log{pair_i}(:,:,cell_i) = logDensityHist;                    
            densHists.norm_lin{pair_i}(:,:,cell_i) = densityHist_norm;
            densHists.norm_log{pair_i}(:,:,cell_i) = logDensityHist_norm;
        end

    end            

        
    [axesLims2D, axesTickLims2D] = deal(zeros(nPairs,4));
    for plot_i = 1:nPairs
        [idx1, idx2] = deal(fetPairs(plot_i,1), fetPairs(plot_i,2));
        x_lims = [fetLims(1, idx1), fetLims(2, idx1)];
        y_lims = [fetLims(1, idx2), fetLims(2, idx2)];
        axesLims2D(plot_i,:) = [x_lims, y_lims];
        x_ticks = x_lims + [1, -1]*diff(x_lims)/(2*nBins);
        y_ticks = y_lims + [1, -1]*diff(y_lims)/(2*nBins);
        axesTickLims2D(plot_i,:) = [x_ticks, y_ticks];
    end
    

end