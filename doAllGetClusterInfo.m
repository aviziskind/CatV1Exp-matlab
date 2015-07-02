function doAllGetClusterInfo

    gids = getAllGids('fd');
    
    for i = 1:length(gids)
        fn = getFileName('prunedClustersStats', gids(i));
        if ~exist(fn, 'file')
        tic; fprintf('%d/%d. Gid = %d.  ', i, length(gids), gids(i)); 
    %         getClusterData(gids(i), 1);        
            exploreTetrodeData(gids(i), 1);        
        toc;
        end
    end


end