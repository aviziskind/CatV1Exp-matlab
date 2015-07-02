function dbUpdateMoviePathNames
    movieTable = 'TBL_MOVIES';
    moviePresTable = 'TBL_MOVIE_PRES';
    hnd = dbOpenExpDb;
    movieIds = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', movieTable);
    
    % update movie name in TBL_MOVIES table
    for mi = 1:length(movieIds)        
        filename = getFieldsFromDatabaseTable(hnd, 'TXT_MOVIE_FILE_NAME', movieTable, {'MOVIE_ID', movieIds(mi)}, [], 1, 1);
        basePath = getName('movieBasePath');
        subPath = getMovieFileType(filename);
        path = [basePath subPath '\'];
        
        updateValueInDatabaseTable(hnd, path, 'TXT_NETWORK_PATH', movieTable, {'MOVIE_ID', movieIds(mi)});
        updateValueInDatabaseTable(hnd, [path filename], 'TXT_MOVIE_FILE', moviePresTable, {'MOVIE_ID', movieIds(mi)});        
                
    end

    

end