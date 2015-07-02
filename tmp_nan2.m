S = load('driftingGratingCells_all.mat');                    
Gids = [S.allCells.Gid];
cellIds = [S.allCells.cellId];                    

% [Gids, m, n] = unique(Gids);
% cellIds = cellIds(m);
for i = 1:length(Gids)
    
    [bins, allHistVals] = dbGetCellSpkStimHists(Gids(i), cellIds(i));
    if any(isnan(allHistVals(:)))
    
        [nOri, nSp, nPh, nCyc, nRep] = dbGetUniqueOriSpPh('Gid', Gids(i), 'length');
        if round(nCyc) == nCyc
            % make sure that # nan in nan columsn is all the same
            n = sum(isnan(allHistVals),3);
            [uVals, count] = uniqueCount(n);
            if length(uVals) ~= 2
                fprintf('Gid = %d. cellId = %d. not unique\n', Gids(i), cellIds(i));
                continue;
            end
            uVal = uVals(2);
            if (uVal ~= nRep)
                fprintf('Gid = %d. more than one\n', Gid);
            end
            3;
        end
    end
end
    

