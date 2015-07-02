function checkMovieNBytesNFrames
    hnd = dbOpenExpDb;

%     movieIdsUsed = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIE_PRES');
%     movieIdsUsed = unique(movieIdsUsed);

    showOKmovies = false;

    [movieIds, moviePath, movieName, nrows, ncols, nfrms] = getFieldsFromDatabaseTable(hnd, ...
        {'MOVIE_ID', 'TXT_NETWORK_PATH', 'TXT_MOVIE_FILE_NAME', 'LNG_N_ROWS', 'LNG_N_COLUMNS', 'LNG_N_FRAMES'}, 'TBL_MOVIES');
    
    for i = 1:length(movieIds)
        f = dir(moviePath{i});
        idx = find(strcmpi({f.name}, movieName{i}), 1);
        if ~isempty(idx)
            
            nBytesExpected = nrows(i)*ncols(i)*nfrms(i);
            nBytesInMovieFile = f(idx).bytes;
            if (nBytesExpected ~= nBytesInMovieFile)
                fprintf([outOf(i, length(movieIds)) ' : (!) For movie %s : Movie file is %d bytes, but should be %d bytes\n'], movieName{i}, nBytesInMovieFile, nBytesExpected );
            end

            txtFileName = strrep(movieName{i}, '.RAW', '.TXT');
            if exist([moviePath{i} txtFileName], 'file');
                param_mtx = load([moviePath{i} txtFileName]);
                
                nFramesInMovieFile = nfrms(i);
                nFramesInTextFile = size(param_mtx, 1);
                
                if (nFramesInTextFile ~= nFramesInMovieFile)
                    fprintf([outOf(i, length(movieIds)) ' : For movie %s : Movie file has %d frames, but text file has %d frames\n'], movieName{i}, nFramesInMovieFile, nFramesInTextFile );
                end
            end
            if showOKmovies
                disp([outOf(i, length(movieIds)) ' : ' movieName{i} ' : ok!']);
            end
            
        else
            disp([outOf(i, length(movieIds)) ' [[ ' movieName{i} ' not found ]]']);
            3;
        end
    
    end
    
    
    
end

%{
updateValueInDatabaseTable(hnd, 8192,  'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', 202, 'TXT_MOVIE_FILE_NAME', 'FGCOMP_2X64X64X8192_1.RAW'});
updateValueInDatabaseTable(hnd, 10240, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', 115, 'TXT_MOVIE_FILE_NAME', 'SPN_16X16X10240_1.RAW'});
updateValueInDatabaseTable(hnd, 10240, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', 143, 'TXT_MOVIE_FILE_NAME', 'MSPN_16X16X10240_1.RAW'});
updateValueInDatabaseTable(hnd, 4800, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', [194, 195, 196, 197]});

%}


