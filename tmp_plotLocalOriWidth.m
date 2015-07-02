function tmp_plotLocalOriWidth
    
    [w_local_B_cell, w_local_B_mu] = getWLocal('driftingGratingCells_GLFcuw8_degree_all_B.mat');
    [w_local_absB_cell, w_local_absB_mu] = getWLocal('driftingGratingCells_GLFcuw8_degree_all_absB.mat');
    [w_local_absBB_cell, w_local_absBB_mu] = getWLocal('driftingGratingCells_GLFcuw8_degree_all.mat');
%     w_local_absB = getWLocal('driftingGratingCells_GLFcuw8_degree_all_absB.mat');

    binE = linspace(0, 200, 71);
    figure(1); clf;
    hist2(w_local_B_cell, w_local_B_mu, binE);
    xlabel('local ori width : unconstrained : +B')

    figure(2); clf;
    hist2(w_local_absB_cell, w_local_absB_mu, binE);
    xlabel('local ori width : +|B| ')

    figure(3); clf;
    hist2(w_local_absBB_cell, w_local_absBB_mu, binE);
    xlabel('local ori width : +|B|+B+{min}')
    
        
    figure(22); clf;
    h = plot(w_local_B_mu, w_local_absB_mu, 'ro', w_local_B_cell, w_local_absB_cell, 'g.');
    hold on; fplot(@(x) x, xlim, 'k:'); 
    legend(h([2,1]), {'cell', 'MU'}, 'location', 'NW')
    xlabel('+B'); ylabel(' +|B|');
    ylim([0 max([w_local_absB_mu(:); w_local_absB_cell(:)])])

    figure(23); clf;
    h = plot(w_local_B_mu, w_local_absBB_mu, 'ro', w_local_B_cell, w_local_absBB_cell, 'g.');
    hold on; fplot(@(x) x, xlim, 'k:'); 
    legend(h([2,1]), {'cell', 'MU'}, 'location', 'NW')
    xlabel('+B'); ylabel(' +|B| + B_{min}');
    ylim([0 max([w_local_absBB_mu(:); w_local_absBB_cell(:)])])
    
    3;
    

    

end


function [w_local_cell, w_local_mu] = getWLocal(filename)

    S = load(filename);
    
    stimTypes = {S.allCells.stimType};
    ori_cells = cellfun(@(s) strncmp(s, 'Grating:Orientation', 15), stimTypes);
    
    oriCells = S.allCells(ori_cells);
    oriCellIds = [oriCells.cellId];
    mu_idx = oriCellIds == 0 | oriCellIds == 100;
    oriStats = [oriCells.stats];
    oriTuningStats = [oriStats.tuningStats];
    
    stats_si = [oriTuningStats.oriStats_si];
    w_local_all = [stats_si.w_ori_local];
    
    w_local_mu = (w_local_all( mu_idx));
    w_local_cell = (w_local_all(~mu_idx));
    
end

function hist2(x1, x2, nbins_arg)
    if length(nbins_arg) == 1
        xlims = lims([nonnans(x1(:)); nonnans(x2(:))]);
        binE = linspace(xlims(1), xlims(2), nbins+1);        
    else
        binE = nbins_arg;
    end
    binC = binEdge2cent(binE);
    vals1 = histcnt(x1, binE);
    vals2 = histcnt(x2, binE);
    bar(binC, [vals1(:), vals2(:)], 1, 'stacked');
    colormap(summer);
    legend('cell', 'MU');
end

