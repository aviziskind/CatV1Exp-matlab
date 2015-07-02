function checkFlashGratingTxtFilesAreCorrect

    S = load('cellsGroups_movie_fg');
    movieGroups = S.movieGroups_fg;
    movieFileNames = {movieGroups.movieFiles};
    movieFileNamesRAW = unique( movieFileNames );
%     movieFileNames = sxtructArrayField(movieGroups, 'movieFiles');
%     movieFileNamesRAW = unique( [movieFileNames{:}]' );
    hnd = dbOpenExpDb;
    for mi = 1:length(movieFileNamesRAW)
        movieIds(mi) = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIES', {'TXT_MOVIE_FILE_NAME', movieFileNamesRAW{mi}});  %#ok<AGROW>
    end









end
