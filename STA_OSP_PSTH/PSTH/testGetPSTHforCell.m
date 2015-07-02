function testGetPSTHforCell 

%     S = load('cellsGroups_movie_long');
%     movieGroups_long = S.movieGroups_long;
    S = load('cellsGroups_movie_fg');
    movieGroups_fg = S.movieGroups_fg;

    idx_long = findInStructArray(movieGroups_fg, 'frameLength_ms', [], @(x) x >= 100);
    longFrameGids = [movieGroups_fg( idx_long ).Gid];

    extFrameLen_ms = 120;

    %     isFlashGratingStimulus = true;
    timeWindow0 = [30 60];
    windowProfile0 = [1 1];

    timeWindowSets = [20 50;
                      25 55;
                      30 60; 
                      35 65
                      40 70;
                      45 75;
                      50 80;
                      55 85;
                      60 90];
    
    modes = {'iterate', 'exhaust'};
    modeNum = 1;
    mode = modes{modeNum};
    
    stab_OSP = 1;
    stab_PSTH = 2;
    stabilityTest = stab_OSP;
    
    nCols = 5;
    nPlots = 2;
    nRowsOfPlots = 2;
    if modeNum == 1
        all_ospSubplotInds  = [(nCols*0) + [2:nCols], (nCols*2) + [1:nCols-1] ];
        all_psthSubplotInds = [(nCols*1) + [2:nCols], (nCols*3) + [1:nCols-1] ];
    else
        all_ospSubplotInds  = [(nCols*0) + [2:nCols], (nCols*2) + [1:nCols] ];
        all_psthSubplotInds = [(nCols*1) + [2:nCols], (nCols*3) + [1:nCols] ];
    end
    osp_ind = all_ospSubplotInds(1);
    psth_ind = all_psthSubplotInds(1);

    function h = subplotOSP( ospdata, idx)
        if idx == 1
            subplot(nPlots*nRowsOfPlots, nCols, (nCols*0)+1);
        elseif idx == 100
            subplot(nPlots*nRowsOfPlots, nCols, (nCols*2) + (nCols) );
        else
            subplot(nPlots*nRowsOfPlots, nCols, osp_ind);
            osp_ind = cycle(osp_ind,  all_ospSubplotInds);
        end
        
        h = imagesc( sum(ospdata, 3) );
        set(gca, 'Position', get(gca, 'OuterPosition'), 'xtick', [], 'ytick', []);            
    end
    
    function h = subplotPSTH( bins, vals, timeWindow, idx)
        if idx == 1
            h = subplot(nPlots*nRowsOfPlots, nCols, (nCols*1)+1 );
        elseif idx == 100
            h = subplot(nPlots*nRowsOfPlots, nCols, (nCols*3) + (nCols) );
        else
            h = subplot(nPlots*nRowsOfPlots, nCols, psth_ind);
            psth_ind = cycle(psth_ind,  all_psthSubplotInds);
        end;

        plotThisPSTH( bins, vals, extFrameLen_ms, [], timeWindow);
        set(gca, 'Position', get(gca, 'OuterPosition'));
    end
    
    figbase = 100*(modeNum-1);
    allDone = [];
    glob_i = gcf+1;

    shortFrameIds = findInStructArray(movieGroups_fg, 'frameLength_ms', [], @(x) x < 17);
    allGids = [movieGroups_fg.Gid];
%     gidsToDo = allGids; 
%     gidsToDo = longFrameGids;
    gidsToDo = [movieGroups_fg(shortFrameIds).Gid];
    
    for gi = 4;%:length(gidsToDo)
        Gid = gidsToDo(gi);
        cellIds = movieGroups_fg(findInStructArray(movieGroups_fg, 'Gid', Gid)).spikingCells;
%         Gid     = movieGroups_long(gi).Gid;
%         cellIds = movieGroups_long(gi).spikingCells;

        for cell_i = 2;%:length(cellIds)
            cellId = cellIds(cell_i);
            if ~isempty(findRows( [Gid, cellId], allDone))
                continue;
            end
            figure; 
%             figure(figbase + glob_i); clf;
            set(gcf, 'name', sprintf('Group %d, cell %d, ', Gid,  cellId) );

            % initial OSP & PSTH
            spkTsRelToFrame_ms = getParsedSpikes('timing', Gid, cellId );
            relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId );
                        
%             nFramesEachPres = cellfun(@length, relContrOfFrameToSpike);
%             nTotalFrames = sum(nFramesEachPres);
            osp_ind = all_ospSubplotInds(1);
            psth_ind = all_psthSubplotInds(1);


            [prevOSP] = getOriSpfPhaseProfile_simple(Gid, relContrOfFrameToSpike);
                subplotOSP( prevOSP, 1);

            [prevPSTH_bins, prevPSTH_vals] = calcPSTH( spkTsRelToFrame_ms, extFrameLen_ms);
            
                subplotPSTH(prevPSTH_bins, prevPSTH_vals, [], 1);
                title({['[' num2str(timeWindow0(1)) ' -- ' num2str(timeWindow0(2)) ']' ], ['Initial PSTH & OSP']});
            
            
            if strcmp(mode, 'iterate')
            
                timeWindow = timeWindow0;
                windowProfile = windowProfile0;

                % iterative PSTH calculation:
                maxIterations = 8;
                timeWindowHistory = zeros(maxIterations+1, 2);
                timeWindowHistory(1,:) = timeWindow;

                iter = 1;
                cc = 0;
                ccThreshold = .98;

                while true;

                    [OSP, oris, sps, phs] = getOSPFromTimeWindow(Gid, cellId,   timeWindow, windowProfile);
                        subplotOSP( OSP, osp_ind );

                    [timeWindow, windowProfile,  PSTH_bins, PSTH_vals] = getPSTHwindowFromOSP(Gid, cellId,   OSP, oris, sps, phs,  extFrameLen_ms);
                        subplotPSTH( PSTH_bins, PSTH_vals, timeWindow, psth_ind);                

                    if stabilityTest == stab_OSP            
                        cc = crossCorrelation(prevOSP, OSP);
                        prevOSP = OSP;            
                    elseif stabilityTest == stab_PSTH 
                        cc = crossCorrelation({prevPSTH_bins, prevPSTH_vals}, {PSTH_bins, PSTH_vals});
                        prevPSTH_bins = PSTH_bins;   prevPSTH_vals = PSTH_vals;
                    end
                    title({['[' num2str(timeWindow(1)) ' -- ' num2str(timeWindow(2)) ']' ], ['(' num2str(iter) ')  ' num2str(cc)]});
                    drawnow;

                    timeWindowStr = ['[Time window: ' num2str(timeWindow(1)) ' -- ' num2str(timeWindow(2)) ']. ' ];
                    if cc >= ccThreshold
                        disp([ '[' num2str(Gid), ', ' num2str(cellId) '] Completed successfully in ' num2str(iter) ' iterations. ' timeWindowStr ]);
                        break;
                    end

                    if ~isempty(findRows( timeWindow, timeWindowHistory(1:iter-1,:) ))
                        disp([ '[' num2str(Gid), ', ' num2str(cellId) '] Algorithm has become cyclic (after ' num2str( iter) ') iterations. ' timeWindowStr ' Aborting ...']);
                        break;
                    end
                    timeWindowHistory(iter+1,:) = timeWindow;

                    iter = iter +1;
                    if iter > maxIterations
                        disp([ '[' num2str(Gid), ', ' num2str(cellId) '] Reached max iterations number (' num2str(iter) '). ' timeWindowStr 'Aborting ...']);
                        break;
                    end
                end

                % check the OSP by using *all* frames now:
                relContrOfFrameToSpike  = getParsedSpikes('frame', Gid, cellId, timeWindow, windowProfile);
                OSP = getOriSpfPhaseProfile_simple(Gid, relContrOfFrameToSpike);
                    subplotOSP( OSP, 100);
                    xlabel({['[' num2str(timeWindow(1)) ' -- ' num2str(timeWindow(2)) ']' ], 'Final OSP'});


            elseif strcmp(mode, 'exhaust')
                ylims2 = zeros(size(timeWindowSets,1),1);
                hnds     = zeros(size(ylims2));
                for tw_i = 1:size(timeWindowSets,1)
                    timeWindowA = timeWindowSets(tw_i,:);
                    
                    [OSP, oris, sps, phs] = getOSPFromTimeWindow(Gid, cellId,  timeWindowA, [1 1]);
                        w = findHowWellTunedOSP(OSP);
                        subplotOSP( OSP, osp_ind ); title( sprintf('w = %3.2f', w) );

                    [timeWindowB, windowProfile, PSTH_bins, PSTH_vals] = getPSTHwindowFromOSP(Gid, cellId,   OSP, oris, sps, phs,  extFrameLen_ms);
                        hnds(tw_i) = subplotPSTH( PSTH_bins, PSTH_vals, timeWindowB, psth_ind);
                        drawVerticalLine(timeWindowA, 'Color', 'c');    
                        title({[' initial: [' num2str(timeWindowA(1)) ' -- ' num2str(timeWindowA(2)) ']' ], ...
                               [' final  : [' num2str(timeWindowB(1)) ' -- ' num2str(timeWindowB(2)) ']' ]});
                   drawnow;
                    ylims = ylim;
                    ylims2(tw_i) = ylims(2);
%                     if tw_i == 9
%                         3;
%                     end
                end
                yl_max = max(ylims2);
                set(hnds, 'ylim', [0, yl_max]);
                
            end
                                                           
            glob_i = glob_i + 1;
            allDone = [allDone; Gid, cellId]; %#ok<AGROW>
        end % of cell
    end % of group

end