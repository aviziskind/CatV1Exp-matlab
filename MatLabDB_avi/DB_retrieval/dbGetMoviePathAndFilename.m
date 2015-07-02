function [path, filename] = dbGetMoviePathAndFilename(hnd, Mid)

    [path, filename] = getFieldsFromDatabaseTable(hnd, {'TXT_NETWORK_PATH', 'TXT_MOVIE_FILE_NAME'}, 'TBL_MOVIES', {'MOVIE_ID', Mid}, [], 1, 1);
       
    % if just 1 output variable - put all info into it.
    if (nargout == 1)
        path = [path filename];
    end
    
end