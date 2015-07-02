function dbFixMovieTableFields
    hnd = dbOpenExpDb;
    movieIdsToCorrect = [202,  115,   143,   194:197];
    nFramesCorrect = [8192, 10240, 10240, 4800*ones(1,4)];
    
    %check whether are incorrect:
    nFramesInDB = zeros(size(movieIdsToCorrect));
    for i = 1:length(movieIdsToCorrect)
        nFramesInDB(i) = getFieldsFromDatabaseTable(hnd, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', movieIdsToCorrect(i)});
        if (nFramesInDB(i) ~= nFramesCorrect(i))
            fprintf('MovieID %d : correcting from %d frames --> %d frames\n', movieIdsToCorrect(i), nFramesInDB(i), nFramesCorrect(i) );
            updateValueInDatabaseTable(hnd, nFramesCorrect(i), 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', movieIdsToCorrect(i)});
        end
    end
    if all(nFramesInDB == nFramesCorrect)
        disp('All frames are correct');
    end
    
end


%     updateValueInDatabaseTable(hnd, 8192,  'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', 202; 'TXT_MOVIE_FILE_NAME', 'FGCOMP_2X64X64X8192_1.RAW'});
%     updateValueInDatabaseTable(hnd, 10240, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', 115; 'TXT_MOVIE_FILE_NAME', 'SPN_16X16X10240_1.RAW'});
%     updateValueInDatabaseTable(hnd, 10240, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', 143; 'TXT_MOVIE_FILE_NAME', 'MSPN_16X16X10240_1.RAW'});
%     updateValueInDatabaseTable(hnd, 4800, 'LNG_N_FRAMES', 'TBL_MOVIES', {'MOVIE_ID', [194, 195, 196, 197]});
