% function comparePSTHSfromTop10AndAll

%     load movieCells_PSTHs;
%     load cellGroups_movie_fg;

    longFrameInds  = findInStructArray(movieGroups_fg, 'frameLength_ms', [], @(x) x > 90);
    groupData = movieGroups_fg(longFrameInds);    
    gridSubPlot(4,5, [101], 2);
    for Gid_i = 2:length(groupData)
        
        Gid = groupData(Gid_i).Gid;
        
        cellIds = groupData(Gid_i).cellIds;
        for cell_i = 1:length(cellIds)
            cellId = cellIds(cell_i);
            
%             [PSTH_bins1, PSTH_vals1] = calcPSTHforLongFrameLengthCells(Gid, cellId);
%             gridSubPlot(1);
%             plotThisPSTH(PSTH_bins1, PSTH_vals1);

            [PSTH_bins2, PSTH_vals2] = calcPSTHforShortFrameFlashGratingCells(Gid, cellId);
            gridSubPlot;
            plotThisPSTH(PSTH_bins2, PSTH_vals2)            
            title([ '(' num2str(Gid_i) ') ' num2str(Gid) ', ' num2str(cellId)  ]);

        end
        
    end
    
        
    


% end
