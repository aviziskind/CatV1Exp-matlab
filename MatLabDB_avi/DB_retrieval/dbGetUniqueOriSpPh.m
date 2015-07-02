function [uOri, uSp, uPh, varargout] = dbGetUniqueOriSpPh(idtype, idval, lengthFlag)
    %{
    % for flashed grating movies:    
        [uOri, uSp, uPh, nReps]         = dbGetUniqueOriSpPh('movieId', movieId)
        [uOri, uSp, uPh, nRepsTotal]    = dbGetUniqueOriSpPh('Gid', Gid);
             % nRepsTotal = #times each frame presented 
        [uOri, uSp, uPh, nPerms, nReps] = dbGetUniqueOriSpPh('Gid', Gid);
             % nPerms = # of different permutations of stimulus order. nReps = # of times each permutation was presented.

    % for flashed grating 'grating' stimuli:
        [uOri, uSp, uPh, nRepsTotal]    = dbGetUniqueOriSpPh('Gid', Gid);
             % nRepsTotal = #times each stimulus type was presented (each presentation will have lasted many frames)

    % for drifting grating ('grating') stimuli:
        [uOri, uSp, uPh, nCyclesTotal]    = dbGetUniqueOriSpPh('Gid', Gid);
             % nRepsTotal = #cycles each grating was drifted

        [uOri, uSp, uPh, nCycles, nReps, uTp] = dbGetUniqueOriSpPh('Gid', Gid);
             % nCycles = # of cycles during each stimulus presentation. nReps = # of times each stimulus was presented. uTf_Hz = temporal frequency. 
%}
    global acceptFgComp
    persistent uosp_data;
%     uosp_data = [];
    redo_all = false;
    
    acceptFgComp = true;
    justReturnLengths = exist('lengthFlag', 'var') && ~isempty(lengthFlag);
    
    if isempty(uosp_data)
        uospFileName = [CatV1Path 'MatLabDB_avi' filesep 'uniqueOriSpPhs.mat'];
        if exist(uospFileName, 'file') && ~redo_all
            S = load(uospFileName);
            uosp_data = S.uosp_data;
        else    
            maxGid = 5336;    
            maxMid = 250;    
            uosp_data.movies = cell(1, maxMid);
            uosp_data.groups = cell(1, maxGid);
            hnd = dbOpenExpDb;
            
            allFgMovieIds = getFieldsFromDatabaseTable(hnd, 'MOVIE_ID', 'TBL_MOVIES', {'TXT_NETWORK_PATH', 'C:\ExperimentDB\Stimuli\Movies\Flashed_Gratings\'});
            allMovieGids = dbGetStimulusGids('movie');
            allMovieGids = allMovieGids(flashedOrDrifting(allMovieGids)==1);
            allGratingGids = dbGetStimulusGids('grating');
            allGratingGids = allGratingGids(flashedOrDrifting(allGratingGids)>0);
            
            doMovies = true;
            doGratings = true;
            
            if doMovies
                % MOVIES
                % 1 - first do each individual fg movie.
                for movieId = allFgMovieIds(:)'
                    [uOri, uSp, uPh, nperms] = getDataForIndividualMovie(movieId);
                    uosp_data.movies{movieId} = {uOri, uSp, uPh, nperms};
                end            

                % 2 - next, do each group, which had multiple movies each
                progressBar('init-', length(allMovieGids), 30);
                disp('Processing movies ... ');
                for Gid = allMovieGids'
                    [uOri, uSp, uPh, nPerms, nReps] = getDataForMovieExperiment(Gid);
                    uosp_data.groups{Gid} = {uOri, uSp, uPh, nPerms, nReps};
                    
                    progressBar;
                end
                progressBar('done');
            end
            
            if doGratings
                disp('Processing gratings ... ');
                progressBar('init-', length(allGratingGids), 30);
                for Gid = allGratingGids'
                    progressBar;                                            
                    [uOri, uSp, uPh, nCyclesEachPres, nRepsEachPres, uTp] = getDataForGratingExperiment(Gid);
                    uosp_data.groups{Gid} = {uOri, uSp, uPh, nCyclesEachPres, nRepsEachPres, uTp};                    
                end                        
                progressBar('done');
            end
            
            save(uospFileName, 'uosp_data');
        end
    end


    switch lower(idtype)
        case {'movieid', 'mid'}
            movieId = idval;
            [uOri, uSp, uPh, nReps] = deal(uosp_data.movies{movieId}{:});
            if justReturnLengths
               [uOri, uSp, uPh] = deal(length(uOri), length(uSp), length(uPh));
            end            
            error(nargoutchk(0,4,nargout));
            varargout = {nReps};
        case {'groupid', 'gid'}
            Gid = idval;
            groupData = uosp_data.groups{Gid};
            [uOri, uSp, uPh, nPerms, nReps] = deal(groupData{1:5});
            uTp = 0; if length(groupData) == 6, uTp = groupData{6}; end;
            if justReturnLengths
               [uOri, uSp, uPh] = deal(length(uOri), length(uSp), length(uPh));
            end            
            error(nargoutchk(0,6,nargout));
            switch  nargout
                case 4, varargout = {nPerms*nReps};
                case 5, varargout = {nPerms, nReps};
                case 6, varargout = {nPerms, nReps, uTp};
            end
        otherwise 
            error('invalid type: first input must be "movieId" or "groupId"')
    end

    
    
    
    function [uOri, uSp, uPh, nPerms] = getDataForIndividualMovie(movieId)
        [phase_deg, sp_pix, ori_deg] = getFlashGratingOSPs_fgmovie(hnd, movieId);
        uOri = unique(ori_deg');  
        uSp = unique(sp_pix');  
        uPh = unique(phase_deg');
        nPerms = length(ori_deg)/(length(uOri)*length(uSp)*length(uPh));    
    end

    function [uOri, uSp, uPh, nUniquePerms, nReps] = getDataForMovieExperiment(Gid)

        movieIds = dbLookup('Mid',  'Gid', Gid);

        [uMovieIds, movieCount] = uniqueCount(movieIds);
        uOri = unique([uosp_data.movies{uMovieIds(1)}{1}]);    
        uSp  = unique([uosp_data.movies{uMovieIds(1)}{2}]);    
        uPh  = unique([uosp_data.movies{uMovieIds(1)}{3}]);    
        for mi = 2:length(uMovieIds)
            assert(isequal(uOri, uosp_data.movies{uMovieIds(mi)}{1}));
            assert(isequal(uSp,  uosp_data.movies{uMovieIds(mi)}{2}));
            assert(isequal(uPh,  uosp_data.movies{uMovieIds(mi)}{3}));
        end
        
        nperms_eachMovie = arrayfun(@(mid) uosp_data.movies{mid}{4}, uMovieIds);  
        assert(length(unique(nperms_eachMovie))==1);
        assert(length(unique(movieCount))==1)
        nRepsTot = sum(movieCount)*nperms_eachMovie(1); 

        nUniquePerms = length(uMovieIds)*nperms_eachMovie(1);
        nReps = nRepsTot/nUniquePerms;                
        
    end    
    
    function [uOri, uSp, uPh, nCyclesEachPres, nRepsEachPres, uTf] = getDataForGratingExperiment(Gid)
        
        Did = dbLookup('Did', 'Gid', Gid);
        [allOri_deg, allSp_pix, allPhase_deg, nSustFrames,  allTempPeriod_frm, allTempPh_deg, gratPresId] = getFieldsFromDatabaseTable(hnd, ...
            {'DBL_ORIENTATION_DEGR', 'DBL_SPATIAL_PERIOD_PIX', 'DBL_SPATIAL_PHASE_DEGR', 'LNG_N_SUSTAINED_FRM', 'DBL_TEMP_PERIOD_FRM', 'DBL_TEMP_PHASE_DEGR', 'GRATING_PRES_ID'}, 'TBL_GRATING_PRES', {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
        %                     [] = getFieldsFromDatabaseTable(hnd, {}, 'TBL_GRATING_PRES', {'DATAFILE_ID', Did});
        %                 [allOri_deg, allSp_pix, allPhase_deg, startTick, endTick] = getFieldsFromDatabaseTable(hnd, {'DBL_ORIENTATION_DEGR', 'DBL_SPATIAL_PERIOD_PIX', 'DBL_SPATIAL_PHASE_DEGR', 'LNG_START_TICK', 'LNG_END_TICK'}, 'TBL_GRATING_PRES', {'DATAFILE_ID', Did});
        uOri = unique(allOri_deg);
        uSp = unique(allSp_pix);
        uPh = unique(allPhase_deg);
                
        if allTempPeriod_frm(1) == 50000  % ie. flashed grating
            
            assert(length(unique(allTempPh_deg))==1)
            assert(length(unique(allTempPeriod_frm))==1);
            
            nreps = length(allOri_deg)/(length(uOri)*length(uSp)*length(uPh));            
            nCyclesEachPres = 1;
            nRepsEachPres = nreps;            
            uTf = 0;
            
            %                 frameLengths = dbConvertTimeMeasures(Did, endTick-startTick, 'tick', 'ms');
            %                 frameLength_ms = mean( frameLengths );
            
        else  %% drifting gratings
            % get phases for each frame
            spatTempPhase_deg = getGratingOSPs_grating(hnd, gratPresId(1));
            uPh = unique(spatTempPhase_deg);
            
            uStims = unique([allOri_deg, allSp_pix, allPhase_deg], 'rows');
            nRepsEachPres = length(allOri_deg)/size(uStims,1);
            if nRepsEachPres ~= round(nRepsEachPres)
                %                             fprintf('Gid = %d: not a whole number of repetitions (%.2f). Needs special attention\n', Gid, nRepsEachPres);
                %                             continue;
                %                             error('Not a whole number of repetitions. Needs special attention');
            end
            nCyclesEachPres = unique(nSustFrames./allTempPeriod_frm);
            %                         if length(unique(nCyclesEachPres))>1
            %                             fprintf('Gid = %d: Presentations uneven\n', Gid);
            %                             continue;
            %                             error('Presentations uneven');
            %                         end
            %                         nCyclesEachPres = nCyclesEachPres(1);

            frmLength_sec = 100/(12*1000);
            tempFreq_Hz = 1./(allTempPeriod_frm * frmLength_sec);
            uTf = unique(tempFreq_Hz);

        end
    end


end


    


    
%             allFG_Gids = find(flashedOrDrifting(1:maxGid)==1);                        
%             progressBar('init-', length(allFG_Gids), 30);
%             for Gid = allFG_Gids'
%                 progressBar;
%                 [allOri_deg, allSp_pix, allPhase_deg, uOri, uSp, uPh] = getOriSpPhaseForEachStimFrame(Gid);
%                 stimType = getStimulusTypeForDid(dbLookup('Did',  'Gid', Gid));
%                 if strcmp(stimType, 'movie')
%                     nrep = length(allOri_deg)/(length(uOri)*length(uSp)*length(uPh));
%                 else
%                     nrep = 1;
%                 end                
%                 uosp_data{Gid,1} = uOri;
%                 uosp_data{Gid,2} = uSp;
%                 uosp_data{Gid,3} = uPh;
%                 uosp_data{Gid,4} = nrep;
%             end

%             for mi = 1:length(uMovieIds)
%                 [phase_deg, sp_pix, ori_deg] = getFlashGratingOSPs_fgmovie(hnd, uMovieIds(mi));
%                 uOri = unique(ori_deg);  uSp = unique(sp_pix);  uPh = unique(phase_deg);
%             
%             
%                 [allOri_deg, allSp_pix, allPhase_deg, uOri, uSp, uPh] = getOriSpPhaseForEachStimFrame(Gid, 'Movie', 'Planned');
%                 nrep = length(allOri_deg)/(length(uOri)*length(uSp)*length(uPh));
%                 frameLength_ms = getFrameLength('Gid', Gid);
%                 
%                 uosp_data(Gid,:) = {uOri, uSp, uPh, nrep, frameLength_ms};
%             end


                % eg 1
                % movieID: M1 M2  M1 M2  M1 M2
                % nreps:   2  2   2  2   2  2
                % here, have have 3 repetitions of 4 permutations
                %
                % eg 2
                % movieID: M1 M2 M3 M4 M5 M6
                % nreps:   2  2  2  2  2  2
                % here, have have 1 repetitions of 12 permutations
