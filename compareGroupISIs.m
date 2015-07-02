function compareGroupISIs(Gid)

    groupSpikes = dbGetSpikes(Gid, [], 'ms');
    spikeTimes = groupSpikes(:,1);
%     cellSpkIds = groupSpikes(:,2);

    spikeSortingData = 'spiker';
%     spikeSortingData = 'klustakwik';
    switch spikeSortingData
        case 'spiker',      cellSpkIds = groupSpikes(:,2);
        case 'klustakwik',
            cluFilename = [CatV1Path 'spikeSorting\KlustaKwik\data\group_' num2str(Gid) '_clusters.mat'];
            S = load(cluFilename);
            cellSpkIds = S.clusterIds;
        otherwise, error('unknown type');
    end    
    
    
    [uCellIds, cellSpikeIdx] = uniqueList(cellSpkIds);
    nCells = length(uCellIds);
    
%     singleCellISIs = cell(1,nCells);
%     doubleCellISIs = cell(nCells,nCells);
    cellISIs = cell(nCells,nCells);
    
    getISIs = @(x) diff(sort(x));
    nTot = nCells*(nCells+1)/2;
    progressBar('init-', nTot);
    for cell_i = 1:nCells        
        progressBar;
        spk_i = spikeTimes(cellSpikeIdx{cell_i});
        cellISIs{cell_i,cell_i} = getISIs( spk_i );
        for cell_j = 1 : cell_i-1
            progressBar;
            spk_j = spikeTimes(cellSpikeIdx{cell_j});
            cellISIs{cell_i, cell_j} = getISIs( [spk_i; spk_j] );
        end
    end
    progressBar('done', nTot);
    
    range_ms = 10;
    nbins = 50;
    binEdges = linspace(0, range_ms, nbins+1);
    binCent = binEdge2cent(binEdges);
    figure(13); clf;
    
    progressBar('init-', nTot);
    for cell_i = 1:nCells
        for cell_j = 1:cell_i
            progressBar;
%             sp_idx = sub2ind([nCells nCells], cell_j, );
            h_ax(cell_i,cell_j) = mySubplot(nCells,nCells, nCells-cell_i+1, cell_j);
            allISIs = cellISIs{cell_i, cell_j};
            isi_cnt = histcnt(allISIs, binEdges);
            h = bar(binCent, isi_cnt, 1, 'edgecolor', 'none');
            3;
            cutOff_ms = 0.8; 
            nSpksBeforeCutoff = nnz(allISIs < cutOff_ms);
%             binsBeforeCutoff = find(binCent < 2);
%             nSpksBeforeCutoff = sum(isi_cnt(binsBeforeCutoff));
            pctSpksBeforeCutoff = nSpksBeforeCutoff / length(allISIs) * 100;
            col = iff(nSpksBeforeCutoff > 10 || pctSpksBeforeCutoff >.5, 'r', 'b');
            prefix = iff(cell_i == cell_j, sprintf('[%d]', uCellIds(cell_i)), sprintf('[%d,%d]', uCellIds(cell_j), uCellIds(cell_i)));
            
            s1 = sprintf('%s<%.1fms;', prefix, cutOff_ms);
            s2 = sprintf('%d spks, (%.1f%%)', nSpksBeforeCutoff, pctSpksBeforeCutoff);
%             h_t(cell_i, cell_j) = title({s1, s2}, 'color', col, 'fontsize', 8);            
            ylims = ylim;
            h_t(cell_i, cell_j) = text(range_ms/2, ylims(2), {s1, s2}, 'color', col, 'fontsize', 8, 'hor', 'center', 'vert', 'top');            
            
            set(gca, 'xtick', [], 'ytick', [])
        end
    end
    
    for cell_i = 1:nCells
        for cell_j = 1:cell_i    
%             tit_pos = get(h_t(cell_i, cell_j), 'position');                        
%             set(h_ax(cell_i, cell_j), 'position', get(h_ax(cell_i, cell_j), 'outerposition'));
%             set(h_t(cell_i, cell_j), 'position', tit_pos);
        end
    end
    
    
end

