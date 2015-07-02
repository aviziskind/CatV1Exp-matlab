function fd_type = flashedOrDrifting(Gid, strFlag)
    persistent fd_data;
%     fd_data = [];
    redo = false;

    if isempty(fd_data)
        flashedDriftingFileName = [CatV1Path 'MatLabDB_avi' filesep 'flashedDrifting.mat'];        
        if exist(flashedDriftingFileName, 'file') && ~redo
            S = load(flashedDriftingFileName);
            fd_data = S.fd_data;
        else
            fd_data = calcAllFlashedOrDriftingData;                                                   
            save(flashedDriftingFileName, 'fd_data');
        end
    end
    
    fd_type = fd_data(Gid);
    if exist('strFlag', 'var') && ~isempty(strFlag)
        strs = {'', 'flashed', 'drifting'};
        fd_type = strs{1+fd_type};
    end
    
end


%----------------------------------------------
function fd_data = calcAllFlashedOrDriftingData
    FLASHED = 1; DRIFTING = 2;    

    maxGid = 5336;    
    fd_data = zeros(maxGid,1);
    hnd = dbOpenExpDb;

    movieGids = dbGetStimulusGids('movie', 'all');
    gratingGids = dbGetStimulusGids('grating', 'all');
    progressBar('init', length(movieGids)+length(gratingGids), 30);
    
    movie_is_flashed = false(size(movieGids));
    movieId = cell(size(movieGids));    
    for Gid_i = 1:length(movieGids)
        progressBar;
        Did = dbLookup('Did', 'Gid', movieGids(Gid_i));
        movieId{Gid_i} = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIE_PRES', {'DATAFILE_ID', Did});
%         movieFileName = getFieldsFromDatabaseTable(hnd, 'TXT_MOVIE_FILE_NAME', 'TBL_MOVIES', {'MOVIE_ID', movieId});
%         movieFileType = getMovieFileType(movieFileName);
        movieFileType = getMovieType('Did', Did);                
        movie_is_flashed(Gid_i) = strcmp(movieFileType, 'Flashed_Gratings');
    end            
    fd_data(movieGids(movie_is_flashed)) = FLASHED;

    %grating stimuli
    grating_is_flashed = false(size(gratingGids));
    grating_is_drifting = false(size(gratingGids));
    fg_stimId = 3;  %getFieldsFromDatabaseTable(hnd, 'STIMULUS_TYPE_ID', 'TBL_STIMULUS_TYPES', {'TXT_STIMULUS_TYPE', 'Flashed Grating Batch'});
    dg_stimIds = [1 2 6 7 9];  % % 'Single Grating', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'
    for i = 1:length(gratingGids)                
        progressBar;
        stimTypeId = dbLookup('Stimulus_id',  'Gid', gratingGids(i));
        grating_is_flashed(i) = stimTypeId == fg_stimId;
        grating_is_drifting(i) = any(stimTypeId == dg_stimIds);
    end            
    progressBar('done');
    fd_data(gratingGids(grating_is_flashed))  = FLASHED;
    fd_data(gratingGids(grating_is_drifting)) = DRIFTING;
end