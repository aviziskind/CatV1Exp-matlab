redoHists = 0;
redoOSPs = 1;

if redoHists
    gids = getAllGids('f');

    progressBar('init=', length(gids));
    for gi = 1:length(gids)    
        S_sorting = load(getFileName('cells', gids(gi)));
        allCellIds = S_sorting.uClustIds;

        fprintf('(%d/%d)', gi, length(gids));
        if any(allCellIds == -1)
            fprintf(' [-1]' );
            [bins, allHistVals] = dbGetCellSpkStimHists(gids(gi), -1);
        end
        if any(allCellIds == 0)        
            fprintf('[0]');
            [bins, allHistVals] = dbGetCellSpkStimHists(gids(gi), 0);
        end
        progressBar(gi);    
        fprintf('\n');
    end
    dbGetCellSpkStimHists('save');
end

if redoOSPs
    gids = getAllGids('f');

    progressBar('init=', length(gids));
    for gi = 1:length(gids)    
        S_sorting = load(getFileName('cells', gids(gi)));
        allCellIds = S_sorting.uClustIds;

        fprintf('(%d/%d)', gi, length(gids));
        if any(allCellIds == -1)
            fprintf(' [-1]' );
            calculatePSTH_STAs_OSP_ForOneCell(gids(gi), -1);
                        
        end
        if any(allCellIds == 0)        
            fprintf('[0]');
            calculatePSTH_STAs_OSP_ForOneCell(gids(gi), 0);
        end
        progressBar(gi);    
        fprintf('\n');
    end
    getOspDataForPsthWindow('save');    
    dbGetCellSpkStimHists('save');    
    
end
3;
