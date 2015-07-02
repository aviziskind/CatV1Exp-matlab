% function generateAllSpikeFiles

    %doing: all flashed grating spike files:
    S = load('cellsGroups_movie_fg');
    movieGroups = S.movieGroups_fg;
    
    nCells = sum( cellfun(@length, {movieGroups.cellIds})) ;
    
    progressBar('init=', nCells);
    for Group_i = 1:length(movieGroups)
        GroupId = movieGroups(Gid_i).Gid;
        cellIds = movieGroups(Gid_i).cellIds;
        
        for ci = 1:length(cellIds)
            filename = getName('spikefile', Gid, cellId);
            if ~exist([path filename ext], 'file')
                generateSpikeFile(Gid, cellId)
            end       
            progressBar;
        end
        
    end
    progressBar('done');        

