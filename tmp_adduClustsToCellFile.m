function tmp_adduClustsToCellFile

    Gids = getAllGids;
    
    progressBar('init-', length(Gids), 40)    
    for i = 1:length(Gids)
        Gid = Gids(i);
        
        sortingFile = getFileName('clustersPruned', Gid, 0);
        S = load(sortingFile);         
        
        [S.uClustIds, S.clustCounts] = uniqueCount(S.clusterIds); 
    
        save(sortingFile, '-struct', 'S');
        progressBar;
    end

    3;

end