function testTunednessAlgorithms(allOSPs)

%{ 
    % when recalculating stats:
    load indivCells_movie_fg
    
    recalculateAllStats;
    generateGratingCellsDatafiles    
%}
    
    
    if nargin < 1
        S = load('flashedGratingCells.mat', 'allOSPs');
        allOSPs = S.allOSPs;
    end
    rankVals = [allOSPs.manualRank];
    
%     mode = 'compareWithRanks'; %vs. sort through visually.    

    showAll = false;
    allOSPs = allOSPs(findInStructArray(allOSPs, 'ori', [], @(x) length(x) > 30));    
    nCells = length(allOSPs);    
    
    % get tunedNess ("w") parameter using present function.
    tic;
    ws = zeros(1, nCells);
    params = [2, .5];
    progressBar('init-', nCells, 30);
    for cell_i = 1:nCells
        progressBar(cell_i);
        show = false;
%         if cell_i == 1141
%             show = true; 
%         end
        ws(cell_i) = findHowWellTunedOSP( allOSPs(cell_i),  params,  show );        
    end 
    progressBar('done');
    toc;
    
    
    % prepare to compare ranks to tunedness parameters
    uW = unique(ws);
    isDiscreteW = (length(uW) / length(ws)) < .1;    

    uRanks = unique(rankVals);
    nRanks = length(uRanks);
    
    ws_forRank = cell(1,nRanks);
    if isDiscreteW        
        nBins = length(uW);
        eps = min(diff(uW))/2;
        w_edges = [uW(:)-eps; max(ws)+.1];        
        
    else
        nBins = 10;
        w_range = [min(ws)-.01, max(ws)+.01];        
        dw = diff(w_range)*.01;
        w_edges = linspace(w_range(1)-dw, w_range(2)+dw, nBins+1);
    end
    w_cents = binEdge2cent(w_edges);
    n = zeros(nRanks, nBins);    
    myRanks = [allOSPs.manualRank];
    for ri = 1:nRanks        
        inds = find(myRanks == uRanks(ri));
        ws_forRank{ri} = ws(inds);
        ord_ws_forRank{ri} = ord(ws_forRank{ri}, 'descend'); %#ok<AGROW>
        origIdxForRankWs{ri} = inds;                %#ok<AGROW>
        
        n_i = histcnt(ws_forRank{ri}, w_edges);
%         n_i(end-1) = n_i(end-1)+n_i(end); n(end) = [];
        n(ri,:) = n_i;
    end
    figure(1);
    bar(w_cents, n'); legend(legendarray('r=', uRanks));
    
    figure(2)
    imagesc(w_cents, 1:nRanks, n); axis xy; colorbar;
    set(gca, 'xtick', w_cents, 'ytick', 1:nRanks); 
    xlabel('W'); ylabel('my rank');
%     bar(n', 'stacked'); legend(legendarray('r=', uRanks));    
    3;
    rw1 = ws_forRank{1}; rw1 = rw1(isfinite(rw1));
    rw4 = ws_forRank{4}; rw4 = rw4(isfinite(rw4));
    [xstar, pErr] = getBestDecisionThreshold(rw1, rw4, 1);
    fprintf('x* = %.3f,  with error prob = %.3f%%\n', xstar, pErr*100);
    
    nJ = 5;
    if showAll 
        gridSubPlot(4,4,  [300 3]);
        figure(300);
        for r_i = 1:nRanks %1:length(myRankInds{1})            
            for j = 1:nJ
                subplot(nJ, nRanks, sub2ind([nRanks nJ], r_i, j));                
                idx = origIdxForRankWs{r_i}(ord_ws_forRank{r_i}(j));                                            
                OSP = allOSPs( idx );
                [groupId, cellId] = deal(OSP.GroupId, OSP.cellId);
                imageOSP( OSP, 'mean:ph', 'SO', 'noLabels', 'noTicks');
%                 colorbar;
                title(sprintf( '[%d,%d] r = %d; w = %g', groupId, cellId, myRanks(idx), ws( idx ) ) );                    
            end                        
        end
        suptitle_2('Highest values of w for each rank');

        figure(310);
        for r_i = 1:nRanks %1:length(myRankInds{1})            
            for j = 1:nJ
                subplot(nJ, nRanks, sub2ind([nRanks nJ], r_i, j));                
                idx = origIdxForRankWs{r_i}(ord_ws_forRank{r_i}(end-j+1));
                OSP = allOSPs( idx );
                [groupId, cellId] = deal(OSP.GroupId, OSP.cellId);
                imageOSP( OSP, 'mean:ph', 'SO', 'noLabels', 'noTicks');
%                 colorbar;
                title(sprintf( '[%d,%d] r = %d; w = %g', groupId, cellId, myRanks(idx), ws( idx ) ) );                    
            end                        
        end        
        suptitle_2('Lowest values of w for each rank');
        
    end


end



        %         gridSubPlot(4,4,  [300 3]);
%         for cell_j = 1:nCells %1:length(myRankInds{1})
%             lastOne = gridSubPlot;
%             idx = myRankInds{1}(cell_j);
%             OSP = allOSPs( idx );
%             [groupId, cellId] = deal(OSP.GroupId, OSP.cellId);
%             imageOSP( OSP, 'mean:ph', 'SO', 'noLabels', 'noTicks');
%             
%             title(sprintf( '(%d) %d, %d.  w = %g', cell_j, groupId, cellId, ws( idx )) );
%             
%             if lastOne, break, end;
%         end
%         
%         gridSubPlot(4,4,  [310 3]);
%         for cell_j = 1:length(myRankInds{4}) %nCells:-1:1
%             lastOne = gridSubPlot;
%             idx = myRankInds{4}(cell_j);
%             OSP = allOSPs( idx );
%             [groupId, cellId] = deal(OSP.GroupId, OSP.cellId);
%             imageOSP( OSP, 'mean:ph', 'SO', 'noLabels', 'noTicks');
%             
%             title(sprintf( '(%d) %d, %d.  w = %g', cell_j, groupId, cellId, ws( idx )));
%             
%             if lastOne, break, end;
%         end
