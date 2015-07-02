function checkSustainedFramesWithMovieLengths
    
    % compareSustainedFramesToLengthOfMovies
    %check if there were any movie presentation where NSustainedFrames >
    %NFramesInMovie (ie. movie was repeated multiple times
    
    hnd = dbOpenExpDb;    
    [Did, nSust, movId] = getFieldsFromDatabaseTable(hnd, {'DATAFILE_ID', 'LNG_N_SUSTAINED_FRM', 'MOVIE_ID'}, 'TBL_MOVIE_PRES');
    
    [uMovIds,m,idx] = unique(movId);
    movLengths = getFieldsFromDatabaseTable(hnd, {'LNG_N_FRAMES'}, 'TBL_MOVIES', {'MOVIE_ID', uMovIds});
    
    movPresLengths = movLengths(idx);
    figure(1);
    hist(nSust - movPresLengths)
        
end