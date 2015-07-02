function [movieType, movieSubType, movieFileNames] = getMovieType(idtype, idval)

    hnd = dbOpenExpDb;
    Did = dbLookup('Did', idtype, idval);
    allMovieFileNames = getFieldsFromDatabaseTable(hnd, 'TXT_MOVIE_FILE', 'TBL_MOVIE_PRES', {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
    movieFileNames = uniqueInOrder(allMovieFileNames);
    nMovieFiles = length(movieFileNames);
    inds = cellfun(@(s) s(end), strfind(movieFileNames, '\')); % remove folder name
    for mi = 1:nMovieFiles
        movieFileNames{mi} = movieFileNames{mi}(inds(mi)+1:end);
    end
    getMovieSubTypes = nargout > 1;
    
    allMovieTypes    = cell(1, nMovieFiles);
    allMovieSubTypes = cell(1, nMovieFiles);
    nRepsFG = zeros(1, nMovieFiles);    

    if ~getMovieSubTypes 
        allMovieSubTypes(:) = {''};
    end
    
    for mi = 1:nMovieFiles
        if getMovieSubTypes
            [allMovieTypes{mi}, allMovieSubTypes{mi}, nRepsFG(mi)] = getMovieFileType(movieFileNames{mi});
        else
            allMovieTypes{mi} = getMovieFileType(movieFileNames{mi});
        end        
    end

    isCounterPhaseMovie = strncmp(movieFileNames{1}, 'CPH_', 4);
    cph_txt = iff(isCounterPhaseMovie, 'C', '');    
    
    [movieType, idx] = unique(allMovieTypes); 
    movieSubType = allMovieSubTypes(idx);
    if length(movieType) == 1
        movieType = movieType{1};
        movieSubType = movieSubType{1};
        if strcmp(movieType, 'Flashed_Gratings')
            Gid = dbLookup('Gid', idtype, idval);
            [uori, usp, uph, nperm, nrep]= dbGetUniqueOriSpPh('Gid', Gid);            
            movieSubType = sprintf('%s(%dx%d)%s', movieSubType, nperm, nrep, cph_txt);
%             movieSubType = [movieSubType '(' num2str(sum(nRepsFG)) ')'];
        end        
    elseif length(movieType) == 2
        % are only two groups like this: both have Noise & FlashGrating.
        %  [Gid = 4360; 4361]
        movieType    = [movieType{1}    '|' movieType{2}];
        movieSubType = [movieSubType{1} '|' movieSubType{2}];
    end

end