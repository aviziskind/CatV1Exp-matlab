function movieFileName = mid_getMovieStimFileName(Gid, frameMode)    
    sd = siteDataFor(Gid);
    allMovieIds = sd.stimulusInfo.movieIds;
    stimType = getGratingStimType(sd.stimType);
    nRepsTotal = stimType.nTrials;
        
%     nStim = length(sd.ori_deg)*length(sd.spPeriod_pix)*length(sd.spPh_deg);
%     nRepsPerMovie = sd.stimulusInfo.nFramesPerPres / nStim;
%     length(allMovieIds) * nRepsPerMovie;    
    
    if strcmp(frameMode, 'all') % considers all movies separately, even if identical
        nReps = nRepsTotal;
    else
        nReps = sscanf(frameMode, '%drep');  
        assert(strcmp(frameMode, sprintf('%drep', nReps)));            
        assert(nReps <= nRepsTotal);
    end
                
    uMovieIds = uniqueListInOrder(allMovieIds);
    switch length(uMovieIds)
        case 1, movieId_str = sprintf('%d', uMovieIds);
        case 2, movieId_str = sprintf('%d_%d', uMovieIds);
        case 16, movieId_str = sprintf('%d-%d', uMovieIds(1), uMovieIds(end));              
        otherwise, error('Unknown movie combination');
    end
    
    nRep_str = sprintf('_%drep', nReps);            
    moviePath = getName('compiledMoviePath');    
    movieFileName = [moviePath 'Movies_' movieId_str nRep_str '.raw'];

end


%{
    if strcmp(frameMode, 'all')
        [uMovieIds, uMovieIdx] = uniqueListInOrder(allMovieIds);
        if length(uMovieIds) == 1
            if length(allMovieIds) == 1
                movieId_str = sprintf('%d', uMovieIds(1));
            else
                movieId_str = sprintf('%d[x%d]', uMovieIds(1), length(uMovieIdx{1}));
            end                

        elseif length(uMovieIds) == 2
            if length(allMovieIds) == 2
                movieId_str = sprintf('%d_%d', uMovieIds);
            else
                nRep = length(allMovieIds)/length(uMovieIds);
                movieId_str = sprintf('[%d_%d]x%d', uMovieIds, nRep);                            
            end
        elseif length(uMovieIds) == 16
            assert(all(diff(uMovieIds)==1)); %ie. 16 consecutive #'s                
            movieId_str = sprintf('[%d--%d]', uMovieIds(1), uMovieIds(end));                
        end

        mode_str = '';
    else
        nReps = sscanf(frameMode, '%drep');
        assert(strcmp(frameMode, sprintf('%drep', nReps)));
        assert(nReps <= stimType.nTrials);
        
        uMovieIds = uniqueListInOrder(allMovieIds);                     
        if length(uMovieIds) < 16
            uMovies_str = cellfun(@num2str, num2cell(uMovieIds), 'un', 0);
            movieId_str = cellstr2csslist(uMovies_str, '_');
        else
            assert(all(diff(uMovieIds)==1)); %ie. 16 consecutive #'s                
            movieId_str = sprintf('[%d--%d]', uMovieIds(1), uMovieIds(end));                
        end
        mode_str = ['_' frameMode];
                
    end
%}


%{
    switch frameMode
        case 'uStim',  
            uMovieIds = uniqueListInOrder(allMovieIds);                     
            if length(uMovieIds) < 16
                uMovies_str = cellfun(@num2str, num2cell(uMovieIds), 'un', 0);
                movieId_str = cellstr2csslist(uMovies_str, '_');
            else
                assert(all(diff(uMovieIds)==1)); %ie. 16 consecutive #'s                
                movieId_str = sprintf('[%d--%d]', uMovieIds(1), uMovieIds(end));                
            end
            mode_str = '_uStim';
            
        case 'uMovie', 
            uMovieIds = uniqueListInOrder(allMovieIds);                                      
            if length(uMovieIds) < 16
                uMovies_str = cellfun(@num2str, num2cell(uMovieIds), 'un', 0);
                movieId_str = cellstr2csslist(uMovies_str, '_');                                
            else
                assert(all(diff(uMovieIds)==1)); %ie. 16 consecutive #'s                
                movieId_str = sprintf('[%d--%d]', uMovieIds(1), uMovieIds(end));                
            end
            mode_str = '';      
            
        case 'all',    movieIds = allMovieIds;

    end
%}

%{
   N       Values
 -----    ----------------
   11  :  112               
   25  :  137                                                                   Gid = 1895
   19  :  138                                                                   Gid = 1857
   14  :  144 
    9  :  145 
   30  :  200_201_200_201_200_201_200_201_200_201_200_201_200_201_200_201       Gid = 2597
    1  :  201_200_201_200                                                       Gid = 2995
   34  :  213_214_213_214_213_214_213_214_213_214_213_214_213_214_213_214       Gid = 4254
   44  :  213_214_215_216_217_218_219_220_221_222_223_224_225_226_227_228       Gid = 4362
   21  :  230_231                                                               Gid = 4462
    4  :  232_233 

%}