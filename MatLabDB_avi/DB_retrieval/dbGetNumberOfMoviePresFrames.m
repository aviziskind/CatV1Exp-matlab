function presFrames = dbGetNumberOfMoviePresFrames(idtype, idval)
    Did = dbLookup('Did',  idtype, idval);
    hnd = dbOpenExpDb;
    movieIds = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIE_PRES', {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
    [uMovieIds, m, idx] = unique(movieIds);
    presFrames = zeros(1, length(movieIds));
    for i = 1 : length(uMovieIds)
        presFrames([idx == i]) = getFieldsFromDatabaseTable(hnd, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', uMovieIds(i)});
    end
end