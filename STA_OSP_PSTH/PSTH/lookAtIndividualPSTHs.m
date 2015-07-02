function lookAtIndividualPSTHs(Gid, cellId)
    [PSTH_bins, PSTH_vals] = getIndividualPSTHs(Gid, cellId, 'OS');
    N = size(PSTH_vals, 2);
	stimIndsInOrder = ord(max(PSTH_vals, [], 1), 'descend');        
    PSTH_vals = PSTH_vals(:,stimIndsInOrder);
    
    N = 12;
    gridSubPlot(3,4, 100)
    clf;
    for i = 1:N
        gridSubPlot; 
        plotThisPSTH(PSTH_bins, PSTH_vals(:,i) )    ;
        title(sprintf('Stimulus # %d', i ));
    end

    
    
end