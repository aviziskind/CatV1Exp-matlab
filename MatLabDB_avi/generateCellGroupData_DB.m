function generateCellGroupData_DB(typesToDo)
    tic;
    hnd = dbOpenExpDb;
    global acceptFgComp;

%     typesAvailable = {'movie', 'grating', 'noise', 'mseq'};
    
    % 1. Movie cells
    movieGroups = struct('Gid', [], 'Did', [], 'presIds', [], 'frameLength_ms', [], 'siteOK', [], 'siteOKdetail', [], ...
        'presOK', [], 'nSyncs', [], 'cellIds', [], 'nSpikes', [], 'locationData', [], 'dataFileInfo', [], 'stimulusInfo', [], 'stimType', [], ...
        'contrast', [], 'isSquare', []);

    % 2. Grating cells 
    gratingGroups = struct('Gid', [], 'Did', [], 'presIds', [], 'frameLength_ms', [], 'siteOK', [], 'siteOKdetail', [], ...
        'presOK', [], 'nSyncs', [], 'cellIds', [], 'nSpikes', [], 'locationData', [], 'dataFileInfo', [], 'stimulusInfo', [], 'stimType', [], ...
        'contrast', [], 'isSquare', [], 'ori_deg', [], 'spPeriod_pix', [], 'spPh_deg', [], 'tempPeriod_sec', [], 'tempPh_deg', []);    
    
    % 3. Noise cells 
    noiseGroups = struct('Gid', [], 'Did', [], 'presIds', [], 'frameLength_ms', [], 'siteOK', [], 'siteOKdetail', [], ...
        'presOK', [], 'nSyncs', [], 'cellIds', [], 'nSpikes', [], 'locationData', [], 'dataFileInfo', [], 'stimType', []);

    % 4. M-sequence cells 
    mseqGroups = struct('Gid', [], 'Did', [], 'presIds', [], 'frameLength_ms', [], 'siteOK', [], 'siteOKdetail', [], ...
        'presOK', [], 'nSyncs', [], 'cellIds', [], 'nSpikes', [], 'locationData', [], 'dataFileInfo', [], 'stimType', []);

    
    
%     matchSpiker = curMatchDB;
    matchSpiker = 1;
    onlyDoFlashedGratingMovies = false;
    
%     groupingType = curGroupingType('');
    groupingType = 'cells';        
%     GroupingType = titleCase(groupingType);
            
    if nargin < 1
        typesToDo = {'movie', 'grating'};
%         typesToDo = {'movie'};
%         typesToDo = {'grating'};
%         typesToDo = {'noise'};
    end

    doMovieGroups = any(strcmp(typesToDo, 'movie'));
    doGratingGroups = any(strcmp(typesToDo, 'grating'));
    doNoiseGroups = any(strcmp(typesToDo, 'noise'));
    doMseqGroups = any(strcmp(typesToDo, 'mseq'));
    
%     if redoMode
%         redoField = 'siteOKdetail';
% %         redoCondition = @(s) ~isempty(s);%~strcmp(s, 'ok');
% %         redoCondition = @(s) strcmp(s, 'warning');
%         redoCondition = @(s) ~isempty(strfind(s, 'sync.VS.db'));
%     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frameNumberFields = {'LNG_N_SUSTAINED_FRM', 'LNG_N_DISPLAYED_FRM', 'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM'};
    baseFieldNames = {'DATAFILE_ID'};
    degPerPixFields = {'DBL_DEGREES_PER_PIXEL', 'LNG_BLOCK_HEIGHT_PIX', 'LNG_BLOCK_WIDTH_PIX'};

%     stimTable = [];
%     noisePresTable = [];    
%     gratingPresTable = [];  %#ok<NASGU>
%     moviePresTable = [];    %#ok<NASGU> 
%     mseqPresTable = [];     %#ok<NASGU>
        
%     ghostGids = [134 470 622 801 802];

    
    function [Did, cellIds, nSpikes, presIds, nFramesPerPres, frameLength_ms, degPerPix, ...
            siteOK, siteOKdetail, presOK, nSyncs, locationData, dataFileInfo, DidIdxs1] = getGroupData(Gid, stimType)
        Did = dbLookup('Did',  'Gid', Gid);
        tablePresName = ['TBL_' upper(stimType) '_PRES'];        
        degPerPix = 'DBL_DEGREES_PER_PIXEL';
        pixPerBlock_h = 'LNG_BLOCK_HEIGHT_PIX';
        pixPerBlock_w = 'LNG_BLOCK_HEIGHT_PIX';
        if matchSpiker
            [cellIds, nSpikes] = dbLookupNumSpikes(Gid);            
        else
            [cellIds, spikeIdxs] = getCellSorting(Gid);
            nSpikes = cellfun(@length, spikeIdxs);
        end
        cellIds = double(cellIds);

%         getStimulusFrameSequence(Gid)
        
        DidIdxs = find(stimTable.DATAFILE_ID == Did); 
        presIds       = stimTable.([upper(stimType) '_PRES_ID'])(DidIdxs);
    
        [siteOK, siteOKdetail, presOK] = testSiteValidity(Gid);
        DidIdxs1 = DidIdxs(1);
        if dbDoesFieldExist(hnd, degPerPix, tablePresName)
            degPerPix = stimTable.(degPerPix)(DidIdxs1);
        else
            degPerPix = 1; % for mseq
        end
        if dbDoesFieldExist(hnd, pixPerBlock_h, tablePresName)
            blockHgt = stimTable.(pixPerBlock_h)(DidIdxs1);
            blockWidth = stimTable.(pixPerBlock_w)(DidIdxs1);
            assert(blockHgt == blockWidth);
        else
            blockHgt = 1; % for mseq, grating
        end        
        degPerPix = degPerPix * blockHgt;        
        
        nFramesDisplayed_all = stimTable.LNG_N_DISPLAYED_FRM(DidIdxs);
        nPreBlankFrames_all = stimTable.LNG_N_PRE_BLANK_FRM(DidIdxs);
        nPrePostFrames_all = stimTable.LNG_N_POST_BLANK_FRM(DidIdxs);                
        nFramesPerPres_all = nFramesDisplayed_all - nPreBlankFrames_all - nPrePostFrames_all;
        
        nFramesPerPres = unique(nFramesPerPres_all);
        if length(nFramesPerPres) > 1
            nFramesPerPres = nFramesPerPres_all;
%             warning(sprintf('Uneven frames : Did = %d\n', Did));
        end
        
        
        frameLength_ms = getFrameLength('Did', Did);        
        presIds = sort( presIds(:) )'; % row vector, sorted (are a handful of cases where are not in order)
        nSyncs = dbLookupNumSyncs(Did);    
        
        LocData = dbGetLocationData('Gid', Gid);
        assert(strcmp(LocData.AnimalType, 'Cat'))
        assert(strcmp(LocData.brainStruct, 'V1'))
        locationData = struct('CatId', LocData.AnimalId, 'PenId', LocData.PenetrId, 'LocId', LocData.LocId, 'ElectrodeType', LocData.ElectrodeType, ...
            'HemiId', LocData.hemiId, 'Hemisphere', LocData.hemi, 'AP', LocData.AP, 'ML', LocData.ML, 'AP_ZERO', LocData.AP_ZERO, 'ML_ZERO', LocData.ML_ZERO, ...
            'AP_offset', LocData.AP - LocData.AP_ZERO, 'ML_offset', LocData.ML - LocData.ML_ZERO, 'depth', LocData.depth );
        
        [dataFileName, nChannels, samplingRate, expDateTime_str_db] = getFieldsFromDatabaseTable(hnd, {'TXT_DATAFILE_NAME', 'LNG_N_CHANNELS', 'DBL_SAMPLING_RATE_HZ', 'DTM_CREATED'}, 'TBL_DATA_FILES', {'DATAFILE_ID', Did});
        dateStrFormat_db = 'm/dd/yyyy HH:MM:SS PM';
        dateStrFormat = 'yyyy/mm/dd HH:MM:SS';
        expDateTime_vec = datevec(expDateTime_str_db, dateStrFormat_db);
        expDateTime_str = datestr(expDateTime_vec, dateStrFormat);
        if isempty(dataFileName)
            error('empty')
        end
        
        nChannelFields = 16;
        channelsUsedFieldnames = arrayfun(@(ch_id) ['BLN_CH_' num2str(ch_id, '%02d')], 1:nChannelFields, 'un', 0);        
        [channelsUsed{1:nChannelFields}] = getFieldsFromDatabaseTable(hnd, channelsUsedFieldnames, 'TBL_GROUPS', {'GROUP_ID', Gid});        
        channelIds = find([channelsUsed{:}]);        
        assert(length(channelIds) >= 3);        
        dfsize = getRawDataFileSize(Did);
        
        bytesPerVal = 2;
        recordingDuration_sec = dfsize / (samplingRate * nChannels * bytesPerVal);        
        dataFileInfo = struct('dataFileName', dataFileName, 'dateCreated', expDateTime_str, 'nChannels', nChannels, 'channelIds', channelIds, 'samplingRateHz', samplingRate, 'filesize', dfsize, 'duration_sec', recordingDuration_sec);

%         eccentricity = getGroupEccentricity(Gid);
        
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1. Movie Groups
    if doMovieGroups
        stimType = 'MOVIE';
        disp('Doing Movie groups');
        
        acceptFgComp = true;
                        
        presIdMissedFrmsFields = {[upper(stimType) '_PRES_ID'], 'LNG_N_MISSED_FRM'};
        movieDetailFields = {'MOVIE_ID', 'LNG_ORIGIN_X', 'LNG_ORIGIN_Y', 'LNG_N_ROWS', 'LNG_N_COLUMNS', 'DBL_MAX_INTENSITY', 'DBL_MIN_INTENSITY'};
        
        allFieldNames = [baseFieldNames, frameNumberFields, presIdMissedFrmsFields, degPerPixFields, movieDetailFields];
        tableName = ['TBL_' stimType '_PRES'];
        stimTable = getTableData(tableName, allFieldNames);
        
        fnames = fieldnames(movieGroups);
        s = [fnames, cell(length(fnames),1)]';
        movieGroups_fg = struct(s{:}, 'ori_deg', [], 'spPeriod_pix', [], 'spPh_deg', [], 'tempPeriod_sec', [], 'tempPh_deg', []);

        movieGids = dbGetStimulusGids(stimType);

        if onlyDoFlashedGratingMovies
            fgIdx = flashedOrDrifting(movieGids) == 1;
            movieGids = movieGids(fgIdx);
        end        
        idxs = 1:length(movieGids);
        movieGroups = repmat( movieGroups, length(movieGids), 1);  % preallocate memory
        movieGroups_fg = repmat( movieGroups_fg, length(movieGids), 1);  % preallocate memory

        
        fprintf('Working on %d movie cell groups ...\n', length(movieGids));
        progressBar('init-', length(movieGids), 60);
        fgInd = 1;
        for i = 1:length(movieGids)

            progressBar(i);
            Gid = movieGids(i);
            [Did, cellIds, nSpikes, presIds, nFramesPerPres, frameLength_ms, degPerPix, siteOK, siteOKdetail, presOK, nSyncs, locData, dataFileInfo] = getGroupData(Gid, stimType);
            
            
            [moviesType, moviesSubType, movieFileNames] = getMovieType('Did', Did);
            moviesSubType = iff( ~isempty(moviesSubType), [':' moviesSubType], '');
            
            Did_idxs = find(stimTable.DATAFILE_ID == Did);
            allMovieIds = stimTable.MOVIE_ID(Did_idxs);
            movieIds = allMovieIds;
            
            nrows = unique(stimTable.LNG_N_ROWS(Did_idxs));
            ncols = unique(stimTable.LNG_N_COLUMNS(Did_idxs));
            int_min = unique(stimTable.DBL_MIN_INTENSITY(Did_idxs));
            int_max = unique(stimTable.DBL_MAX_INTENSITY(Did_idxs));            
            origin_x = unique(stimTable.LNG_ORIGIN_X(Did_idxs));
            origin_y = unique(stimTable.LNG_ORIGIN_Y(Did_idxs));
            pixPerBlock = unique(stimTable.LNG_BLOCK_HEIGHT_PIX(Did_idxs));
            
            screenH_pix = 768;
            screenW_pix = 1024;
            assert(length([nrows, ncols, int_min, int_max, origin_x, origin_y])==6);
            
            stimulusInfo = struct('movieIds', movieIds, 'movieFiles', {movieFileNames}, 'nFramesPerPres', nFramesPerPres, ...
                    'origin_x', origin_x, 'origin_y', origin_y, 'pixPerBlock', pixPerBlock, 'nrows', nrows, 'ncols', ncols, ...
                    'degreesPerBlock', degPerPix, 'screenH_pix', screenH_pix, 'screenW_pix', screenW_pix, 'int_min', int_min, 'int_max', int_max );
%             movieIds = unique( allMovieIds );
            
            movieGroup_i = struct('Gid', Gid, 'Did', Did, 'presIds', presIds, 'frameLength_ms', frameLength_ms, ...
                'siteOK', siteOK, 'siteOKdetail', siteOKdetail, 'presOK', presOK, 'nSyncs', nSyncs, 'cellIds', cellIds, 'nSpikes', nSpikes, ...
                'locationData', locData, 'dataFileInfo', dataFileInfo, 'stimulusInfo', stimulusInfo, 'stimType', ['Movie:' moviesType moviesSubType], ...
                'contrast', 100, 'isSquare', nan);                        
  %         isSquare', [], 'ori_deg', [], 'spPeriod_pix', [], 'spPh_deg', [], 'tempPeriod_sec', [], 'tempPh_deg', []);            
            
            if strcmp(moviesType, 'Flashed_Gratings')                
                movieGroup_i.isSquare = false;
                
                s = [fnames, struct2cell(movieGroup_i)]';
                [uori, usp, uph] = dbGetUniqueOriSpPh('Gid', Gid);
                movieGroups_fg(fgInd) = struct(s{:}, 'ori_deg', uori, 'spPeriod_pix', usp, 'spPh_deg', uph, 'tempPeriod_sec', inf, 'tempPh_deg', 0);
                fgInd = fgInd+1;
            end                        
            movieGroups(idxs(i)) = movieGroup_i;
            
        end
        progressBar('done');        
        
        movieGroups_all = movieGroups;        
        movieGroups_fg_all = movieGroups_fg;
        % all movies:        
        if ~onlyDoFlashedGratingMovies
            
            idx_withSpikes = arrayfun(@(s) ~isempty(s.cellIds), movieGroups_all);
            movieGroups_all = movieGroups_all(idx_withSpikes); % the ones without spikes are just duplicates.        
            
            idx_withSiteOK = arrayfun(@(s) strcmp(s.siteOK, 'ok'), movieGroups_all);
            movieGroups = movieGroups_all(idx_withSiteOK);         %#ok<NASGU>
        end
        
        % flashed grating movies:
        fgIdx_withSpikes = arrayfun(@(s) ~isempty(s.cellIds), movieGroups_fg_all);        
        movieGroups_fg_all = movieGroups_fg_all(fgIdx_withSpikes);
        % there are 4 sites with slightly uneven frames (Gid = 3073, 4768, 5100, 5146), and this
        % generates a warning, but the effects on the analysis are so negligible that they can still
        % be included in the analysis.
        fgIdx_withSiteOK = arrayfun(@(s) strcmp(s.siteOK, 'ok') || strcmp(s.siteOK, 'warning'), movieGroups_fg_all);        
        movieGroups_fg = movieGroups_fg_all( fgIdx_withSiteOK );  %#ok<NASGU>
        
        if ~onlyDoFlashedGratingMovies                   
            save([CatV1Path groupingType 'Groups_movie_all_DB.mat'], 'movieGroups_all', 'movieGroups_fg_all');
            save([CatV1Path groupingType 'Groups_movie_DB.mat'],     'movieGroups');
        end
        save([CatV1Path groupingType 'Groups_movie_fg_all_DB.mat'],  'movieGroups_fg_all');
        save([CatV1Path groupingType 'Groups_movie_fg_DB.mat'],  'movieGroups_fg');
        disp('Completed Movie Groups');
    end
    
    % 2. Grating Groups
    if doGratingGroups
        stimType = 'GRATING';
        disp('Doing Grating groups');        

        gratingGids = dbGetStimulusGids(stimType);
        idxs = 1:length(gratingGids);
        gratingGroups = repmat( gratingGroups, length(gratingGids), 1);  % preallocate memory
        
        gratingParamFields = {'DBL_ORIENTATION_DEGR', 'DBL_TEMP_PHASE_DEGR', 'DBL_TEMP_PERIOD_FRM', 'DBL_SPATIAL_PHASE_DEGR', 'DBL_SPATIAL_PERIOD_PIX', 'BLN_SQUARE_GRATING', 'DBL_MIN_INTENSITY', 'DBL_MAX_INTENSITY'};
        presIdMissedFrmsFields = {[upper(stimType) '_PRES_ID'], 'LNG_N_MISSED_FRM'};
        allFieldNames = [baseFieldNames, frameNumberFields, gratingParamFields, presIdMissedFrmsFields, degPerPixFields(1)];
        tableName = ['TBL_' stimType '_PRES'];
        stimTable = getTableData(tableName, allFieldNames);
        
        all_nCycles = zeros(1,length(gratingGids));
        all_nReps = zeros(1,length(gratingGids));
        
        fprintf('Working on %d grating cell groups ...\n', length(gratingGids));
        progressBar('init-', length(gratingGids), 60);
        for i = 1:length(gratingGids)
            progressBar(i); 
            Gid = gratingGids(i);
            [Did, cellIds, nSpikes, presIds, nFramesPerPres, frameLength_ms, degPerPix, siteOK, siteOKdetail, ...
                presOK, nSyncs, locData, dataFileInfo, DidIdx1] = getGroupData(Gid, stimType);
            
            isSquareVal = stimTable.('BLN_SQUARE_GRATING')(DidIdx1);
            maxIntensity = stimTable.('DBL_MAX_INTENSITY')(DidIdx1);
            
            contrastPct = round( (maxIntensity-(1-maxIntensity))*100 );
            isSquare = isSquareVal~=0;
            sqrStr = iff(isSquare, ']', ')');
            [stimulusType, gratingType] = getStimulusTypeForDid(Did);
            [uOri, uSp, uPh, nCycles, nRep] = dbGetUniqueOriSpPh('Gid', Gid);

            all_nCycles(i) = nCycles(1);
            all_nReps(i) = nRep;            
            
            makeStr = @(x) iff(length(x) ==1, num2str(x), ['[' cellstr2csslist(arrayfun(@num2str, x, 'un', 0)) ']']);
            nStimStr = sprintf('%dx%dx%d', length(uOri), length(uSp), length(uPh));
            nCyclesStr = makeStr(nCycles);
            nRepStr = makeStr(nRep);
                
            stimTypeStr = [stimulusType ':' gratingType ':' nStimStr '(' nCyclesStr 'x' nRepStr ')' ':' sqrStr ];
                % gratingTypes = {'Single Grating', 'Orientation Batch', 'Flashed Grating Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'};
            
            DidIdxs = find(stimTable.DATAFILE_ID == Did);
%             ori_deg = stimTable.('DBL_ORIENTATION_DEGR')(DidIdxs);
            tempPh_deg = stimTable.('DBL_TEMP_PHASE_DEGR')(DidIdxs);
            tempPeriod_frm = stimTable.('DBL_TEMP_PERIOD_FRM')(DidIdxs);
            int_min = unique(stimTable.DBL_MIN_INTENSITY(DidIdxs));
            int_max = unique(stimTable.DBL_MAX_INTENSITY(DidIdxs));            

%             spPh_deg = stimTable.('DBL_SPATIAL_PHASE_DEGR')(DidIdxs);
%             spPeriod_pix = stimTable.('DBL_SPATIAL_PERIOD_PIX')(DidIdxs);
            
%             ori_deg = unique(ori_deg');
%             tempPh_deg = unique(tempPh_deg');
            tempPeriod_frm = unique(tempPeriod_frm');
            if (tempPeriod_frm >= 50000), tempPeriod_frm = inf; end
            tempPeriod_sec = dbConvertTimeMeasures(Did, tempPeriod_frm, 'frame', 'sec');
            if ~isinf(tempPeriod_sec)
                3;
            end
%             spPh_deg = unique(spPh_deg');
%             spPeriod_pix = unique(spPeriod_pix');            

            stimulusInfo = struct('nFramesPerPres', nFramesPerPres, ...
                     'degreesPerBlock', degPerPix, 'int_min', int_min, 'int_max', int_max );
            
            
            gratingGroups(idxs(i)) = struct('Gid', Gid, 'Did', Did, 'presIds', presIds, 'frameLength_ms', frameLength_ms, ...
                'siteOK', siteOK, 'siteOKdetail', siteOKdetail, 'presOK', presOK, 'nSyncs', nSyncs, 'cellIds', cellIds, 'nSpikes', nSpikes, ...
                'locationData', locData, 'dataFileInfo', dataFileInfo, 'stimulusInfo', stimulusInfo, 'stimType', stimTypeStr,  ...
                'contrast', contrastPct, 'isSquare', isSquare, 'ori_deg', uOri(:)', 'spPeriod_pix', uSp(:)', 'spPh_deg', uPh(:)', 'tempPeriod_sec', tempPeriod_sec, 'tempPh_deg', tempPh_deg);
        end 
  %         isSquare', [], 'ori_deg', [], 'spPeriod_pix', [], 'spPh_deg', [], 'tempPeriod_sec', [], 'tempPh_deg', []);                    
        
        gratingGroups_all = gratingGroups; 
        idxWithSpikes = arrayfun(@(s) ~isempty(s.cellIds), gratingGroups_all); % the ones without spikes are just duplicates.        
        gratingGroups_all = gratingGroups_all(idxWithSpikes); 
        gratingGroups = gratingGroups_all;  % _all contains everything except the ones just discarded (duplicates with no spikes)
        
        % remove the handful with non-integer number of cycles or repetitions
        wholeNumberCyclesReps = ((all_nCycles == round(all_nCycles)) & (all_nReps == round(all_nReps)));
        wholeNumberCyclesReps = wholeNumberCyclesReps(idxWithSpikes);                         
        gratingGroups = gratingGroups(wholeNumberCyclesReps);                        

        % remove ones where were problems with presentation
        idxWithSiteOK = arrayfun(@(s) strcmp(s.siteOK, 'ok'), gratingGroups);        
        gratingGroups = gratingGroups(idxWithSiteOK);        

        % remove square gratings;
        idx_square = ([gratingGroups.isSquare] == 1);
        gratingGroups(idx_square) = [];         
        
        % only keep gratings with contrast at 100% (this is 95% of them anyway)
        idx_contrast100 = ([gratingGroups.contrast] == 100);
        gratingGroups = gratingGroups(idx_contrast100);
        
        allGratStimTypes = {gratingGroups.stimType};        
        % only keep for dOr and dSf, drifting gratings at 2 or 3 Hz (this is 95% of them anyway)
        tf_Hz = {gratingGroups.tempPeriod_sec};
        ori_spat_batch = cellfun(@(s) any(strncmp(s, {'Grating:Ori', 'Grating:Spa'}, 11)), allGratStimTypes);        
        idx_2_3_Hz = cellfun(@(tF) (1/tF(1) == 2 | 1/tF(1) == 3), tf_Hz);
        idx_keep = ~ori_spat_batch | (ori_spat_batch & idx_2_3_Hz);  % only discard non-2/3 Hz if in 'orientation batch' or 'spat freq batch'.

        gratingGroups = gratingGroups(idx_keep);        
        allGratStimTypes = {gratingGroups.stimType};        
        
        fOrSphInds = cellfun(@(s) strncmp(s, 'Grating:Flashed Grating Batch', 15), allGratStimTypes);
        dSngInds   = cellfun(@(s) strncmp(s, 'Grating:Single Grating', 15), allGratStimTypes);
        dOrInds    = cellfun(@(s) strncmp(s, 'Grating:Orientation Batch', 15), allGratStimTypes);
        dSfInds    = cellfun(@(s) strncmp(s, 'Grating:Spatial Frequency Batch', 15), allGratStimTypes);
        dTfInds    = cellfun(@(s) strncmp(s, 'Grating:Temporal Frequency Batch', 15), allGratStimTypes);
        dFreeInds  = cellfun(@(s) strncmp(s, 'Grating:Free Grating Batch', 15), allGratStimTypes);
                        
        gratingGroups_fOrSph = gratingGroups( fOrSphInds ); %#ok<NASGU> 
        gratingGroups_dSng   = gratingGroups( dSngInds   ); %#ok<NASGU> 
        gratingGroups_dOr    = gratingGroups( dOrInds    ); %#ok<NASGU> 
        gratingGroups_dSf    = gratingGroups( dSfInds    ); %#ok<NASGU>
        gratingGroups_dTf    = gratingGroups( dTfInds    ); %#ok<NASGU>        
        gratingGroups_dFree  = gratingGroups( dFreeInds  ); %#ok<NASGU>        
                 
        save([CatV1Path groupingType 'Groups_grating_all_DB.mat'],   'gratingGroups_all');
        save([CatV1Path groupingType 'Groups_grating_DB.mat'],       'gratingGroups');
        save([CatV1Path groupingType 'Groups_grating_fOrSph_DB.mat'],'gratingGroups_fOrSph'); % 'Flashed Grating' (varies in Ori/Sph)
        save([CatV1Path groupingType 'Groups_grating_dSng_DB.mat'],  'gratingGroups_dSng');   % 'Single grating' (single ori/sp)
        save([CatV1Path groupingType 'Groups_grating_dOr_DB.mat'],   'gratingGroups_dOr');    % 'Orientation Batch' (1 spf, varies in orientation)
        save([CatV1Path groupingType 'Groups_grating_dSf_DB.mat'],   'gratingGroups_dSf');    % 'Spatial Frequency Batch' (a few oris, many SpatFreq)
        save([CatV1Path groupingType 'Groups_grating_dTf_DB.mat'],   'gratingGroups_dTf');    % 'Temporal Frequency Batch' (varies in TempFreq)
        save([CatV1Path groupingType 'Groups_grating_dFree_DB.mat'], 'gratingGroups_dFree');  % 'Free Grating Batch' (varies in Ori/Spf/TempFreq)

        disp('Completed Grating Groups');
    end
    
    
    % 3. Noise Groups
    if doNoiseGroups
        stimType = 'NOISE';
        disp('Doing noise groups');
        
        presIdMissedFrmsFields = {[upper(stimType) '_PRES_ID'], 'LNG_N_MISSED_FRM'};
        allFieldNames = [baseFieldNames, frameNumberFields, presIdMissedFrmsFields, degPerPixFields];
        tableName = ['TBL_' stimType '_PRES'];
        stimTable = getTableData(tableName, allFieldNames);
        
        noiseGids = dbGetStimulusGids(stimType);
        idxs = 1:length(noiseGids);
        noiseGroups = repmat( noiseGroups, length(noiseGids), 1);  % preallocate memory
        
        fprintf('Working on %d noise %s groups ...\n', length(noiseGids), groupingType);
        progressBar('init-', length(noiseGids), 60);
        for i = 1:length(noiseGids)
            progressBar(i);
            Gid = noiseGids(i);            
            [Did, cellIds, nSpikes, presIds, nFramesPerPres, frameLength_ms, degPerPix, siteOK, siteOKdetail, presOK, nSyncs, locData, dataFileInfo] = getGroupData(Gid, stimType);

            noiseGroups(idxs(i)) = struct('Gid', Gid, 'Did', Did, 'presIds', presIds, 'frameLength_ms', frameLength_ms, 'degPerPix', degPerPix, ...
                'siteOK', siteOK, 'siteOKdetail', siteOKdetail, 'presOK', presOK, 'nSyncs', nSyncs, 'cellIds', cellIds, 'nSpikes', nSpikes, ...
                'locationData', locData, 'dataFileInfo', dataFileInfo, 'stimType', 'Noise');
        end       
        noiseGroups_all = noiseGroups;
        
        idx_withSpikes = arrayfun(@(s) ~isempty(s.cellIds), noiseGroups_all); 
        idx_withSiteOK = arrayfun(@(s) strcmp(s.siteOK, 'ok'), noiseGroups_all);
        noiseGroups = noiseGroups_all(idx_withSpikes & idx_withSiteOK); %#ok<NASGU> 
                        
        save([CatV1Path groupingType 'Groups_noise_all_DB.mat'], 'noiseGroups_all');
        save([CatV1Path groupingType 'Groups_noise_DB.mat'],     'noiseGroups');
        
        disp('Completed Noise Groups');
    end


    % 4. M-sequence Groups
    if doMseqGroups
        stimType = 'MSEQ';
        disp('Doing M-sequence groups');
        mseqGids = dbGetStimulusGids(stimType);

        presIdMissedFrmsFields = {[upper(stimType) '_PRES_ID'], 'LNG_N_MISSED_FRM'};
        allFieldNames = [baseFieldNames, frameNumberFields, presIdMissedFrmsFields];
        tableName = ['TBL_' stimType '_PRES'];
        stimTable = getTableData(tableName, allFieldNames);
        
        mseqGroups = repmat( mseqGroups, length(mseqGids), 1); % preallocate memory
                
        progressBar('init-', length(mseqGids), 60);
        for i = 1:length(mseqGids)
            progressBar(i);
            Gid = mseqGids(i);
            [Did, cellIds, nSpikes, presIds, nFramesPerPres, frameLength_ms, degPerPix, siteOK, siteOKdetail, presOK, nSyncs, locData, dataFileInfo] = getGroupData(Gid, stimType);
            
            mseqGroups(i) = struct('Gid', Gid, 'Did', Did, 'presIds', presIds, 'frameLength_ms', frameLength_ms, 'degPerPix', degPerPix, ...
                'siteOK', siteOK, 'siteOKdetail', siteOKdetail, 'presOK', presOK, 'nSyncs', nSyncs, 'cellIds', cellIds, 'nSpikes', nSpikes, ...
                'locationData', locData, 'dataFileInfo', dataFileInfo, 'stimType', 'Mseq');
        end        
        mseqGroups_all = mseqGroups; 
                
        idx_withSpikes = arrayfun(@(s) ~isempty(s.cellIds), mseqGroups_all); 
        idx_withSiteOK = arrayfun(@(s) strcmp(s.siteOK, 'ok'), mseqGroups_all);
        mseqGroups = mseqGroups_all(idx_withSpikes & idx_withSiteOK); %#ok<NASGU>         
        
        save([CatV1Path groupingType 'Groups_mseq_all_DB.mat'], 'mseqGroups_all');
        save([CatV1Path groupingType 'Groups_mseq_DB.mat'],     'mseqGroups');
        disp('Completed Msequence Groups');
    end

    
    toc;
end


%         if all( dbDoesFieldExist(hnd, frameNumberFields, tableName ) );
%             
%             nSustainedFrames = stimTable.LNG_N_SUSTAINED_FRM(DidIdxs);
%             nDisplayedFrames = stimTable.LNG_N_DISPLAYED_FRM(DidIdxs);
%             preBlankFrames   = stimTable.LNG_N_PRE_BLANK_FRM(DidIdxs);
%             postBlankFrames  = stimTable.LNG_N_POST_BLANK_FRM(DidIdxs);
%                         
%             presOK = [(nMissedFrames == 0) & (nDisplayedFrames == nSustainedFrames + preBlankFrames + postBlankFrames) ]';
%         else
%             presOK = (nMissedFrames == 0)';
%         end



%         % compile flashed grating movies, with extra fg-specific data:
%         fprintf('Compiling flashed grating groups... '); tic;
%         fg_inds = findInStructArray(movieGroups_all, 'stimType', [], @(s) strncmp(s, 'Movie:Flashed_Grating', 19));
%         for i = 1:length(fg_inds)
%             s = [fnames, struct2cell(movieGroups_all(fg_inds(i)))]';
%             indMovieFileNames = strcmp(fnames, 'movieFiles');
%             indGid = strcmp(fnames, 'Gid');
%             Gid = s{2, indGid};
%             movieFileNames = s(2, indMovieFileNames);
%             s(:,indMovieFileNames) = [];
% %             s = s';
%             [allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph] = getOriSpPhaseForEachStimFrame(Gid, 'Movie', 'Planned');
%             movieGroups_fg_all(i) = struct(s{:}, 'movieFiles', movieFileNames, 'ori_deg', uori, 'spPeriod_pix', usp, 'spPh_deg', uph, 'tempPeriod_sec', inf, 'tempPh_deg', 0);            
%         end
%         fprintf('done \n'); toc;
%         indsWithSiteOK = findInStructArray(movieGroups_fg_all, 'siteOK', [], @(s) strcmp(s, 'ok') );                
%         movieGroups_fg = movieGroups_fg_all(indsWithSiteOK);


%         vname_all = ['grating' GroupingType 'Groups_all'];
%         vname     = ['grating' GroupingType 'Groups'];
%         vname_fOrSph = ['grating' GroupingType 'Groups_fOrSph'];
%         vname_dSng     = ['grating' GroupingType 'Groups_dSng'];
%         vname_dOr     = ['grating' GroupingType 'Groups_dOr'];
%         vname_Sf     = ['grating' GroupingType 'Groups_dSf'];
%         


%         movieGids = [ 1194 1195 1201 1207 1215 1432 1433 1441 1442 1444 1474 1491 1493 1495 1499 1515 1521 1531 1534 1535 ...
%                         1536 1662 1709 1790 1861 1983 2132 2138 2139 2158 2159 2162 2165 2166 2170 2178 2238 2250 2460 2601 ...
%                         2625 2631 2655 2657 2665 2677 2715 2723 2737 2747 2775 2791 2805 2855 2887 2921 2957 2981 3035 3047 ...
%                         3051 3059 3061 3073 3077 3091 3093 3125 3151 3215 3373 3397 3415 3467 3637 3709 3719 3729 3731 3775 ...
%                         3841 3997 4023 4051 4117 4127 4312 4344 4362 4768 5100 5146 ];
