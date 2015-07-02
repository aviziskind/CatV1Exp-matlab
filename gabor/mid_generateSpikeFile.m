function [spksPerFrame_exact, spksPerFrame_uint8] = mid_generateSpikeFile(Gid, cellId, windowMode, trialMode, frameMode, responseType, varargin) 
%     mid_generateSpikeFile(Gid, cellId, timeWindow_ms, windowProfile_ms, frameMode, timeWindow_mode);
%     mid_generateSpikeFile(Gid, cellId, frameMode, trialMode) 
%     mid_generateSpikeFile(Gid, cellId, frameMode, trialMode, relContrOfFrameToSpike) 
%     mid_generateSpikeFile(Gid, cellId, frameMode, trialMode, spkWindow_ms, spkProfile) 
%     mid_generateSpikeFile('all') 
    

    if strcmp(Gid, 'all')
%         S = load('flashedGratingCells_GLFcuw8_phase.mat');
%         gids = [S.allCells.Gid];
%         cellids = [S.allCells.cellId];
%         idx = cellids > 0; gids = gids(idx); cellids = cellids(idx);        
        [allGids, allCellIds] = getAllGids('f');        
%         allGids = 4470;
%         allCellIds = 2;
        
        responseType = curResponseType('');
        doAllFrameModes = 0;
%         allGids = 4502;
%         allCellIds = 4;
        nCells = length(allCellIds);        
        
        useDataFromIndivFile = 0;
        
        S2 = load('indivCells_GLFcuw8_movie_fg.mat');  
        progressBar('init-', nCells, 60)
        for ci = 1:nCells
            Gid = allGids(ci);
            if Gid == 4482
                3;
            end
            cellId = allCellIds(ci);
            nm = getName('celldata', Gid, cellId);
            v = S2.(nm);
            assert(v.Gid == Gid && v.cellId == cellId);

            fprintf(' == %d/%d: Gid = %d, cellId = %d == \n', ci, nCells, Gid, cellId);
            
            trialModes = {'all', 'odd', 'even'};            
            
%             windowModes = {[29+1/6, 62.5], [58+1/3, 91+2/3]};
% %             allTimeWindows = {'best', [29, 62], [58, 91]};
            
% 			windowModes = {[30 60], [60 90]};         
            windowModes = {[29, 62], [58, 91]};         
            
% 			windowModes = {'best'};            
			doExtraWindows = 1; %isfield(v.STAs, 'fixedWindowSTAs');
% 			if doExtraWindows
% 				extraWindows_ms = v.STAs.fixedWindowSTAs.windows_ms;
% 				nExtra = size(extraWindows_ms, 1);
% 				extraWindows = mat2cell( extraWindows_ms, ones(nExtra, 1) , 2)';
% 				windowModes = [windowModes, extraWindows]; %#ok<AGROW>
% 			end
			
			nWindows = length(windowModes);
			for wind_i = 1:nWindows
                windowMode_i = windowModes{wind_i};
                curTimeWindow(windowMode_i);
                PSTHdata = getPSTHforCell(Gid, cellId);
                if ischar(windowMode_i) && strcmp(windowMode_i, 'best')
%                     if useDataFromIndivFile
%                         PSTHdata = v.PSTH;
%                     else
%                         PSTHdata = getPSTHforCell(Gid, cellId);
%                     end

%                     spkWindow_ms = psth.timeWindow_ms;
                    [L_bin, R_bin] = dealV( PSTHdata.timeWindow_bins );

                    windowProfile = PSTHdata.windowProfile;
                    spkWindow_ms = PSTHdata.timeWindow_ms;
                else
                    %%
                    [L_bin, R_bin] = getLRbins_windowprofile(windowModes{wind_i}, PSTHdata);
                    
% 					spkWindow_ms_approx = windowModes{wind_i};
%                     
%                     dBin = (1000/120)/2; binCents = [-300 + dBin/2 : dBin : 200];
%                     L_bin = indmin(abs( binCents - spkWindow_ms_approx(1))); 
%                     R_bin = indmin(abs( binCents - spkWindow_ms_approx(2))); 
%                                        
%                     spkWindow_ms = binCents([L_bin, R_bin]) + dBin/2*[-1, 1];
					windowProfile = [];
                end
%                 if strcmp(responseType, 'raw')  || true;
%     %             relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, spkWindow_ms, windowProfile );
%                     relContrOfFrameToSpike = getSpkResponseToEachFrame(Gid, cellId, spkWindow_ms, windowProfile, 'stimOrder');
%                     responseToEachStim = relContrOfFrameToSpike;
%                 end
                
                if strcmp(responseType, 'gainCorrected')  || true;
                
                    osp_full = double( getOspDataForPsthWindow(Gid, cellId, [], [], L_bin, R_bin, windowProfile, 'osp_full') );
                    [nOri, nSpf, nPh, nTrials] = size(osp_full); 
                    responseToEachStim = reshape(osp_full, [nOri * nSpf * nPh, nTrials])';
%                     alldiffs = responseToEachStim2(:)/nanmean(responseToEachStim2(:)) - responseToEachStim(:)/nanmean(responseToEachStim(:));
%                     assert(max(abs(alldiffs)) < 1e-4);
                    
                end
                
                if doAllFrameModes
                    stimType = getGratingStimType(Gid);
                    nTrials = stimType.nTrials;
                    if nTrials == 16
                        frameModes = {'1rep', '2rep', '4rep', 'all'};
                    elseif nTrials == 10
                        frameModes = {'1rep', '2rep', '5rep', 'all'};
                    else 
                        frameModes = {'1rep', '2rep', 'all'};
                    end
                else
                    frameModes = {'1rep'};
                end
				nNrep = length(frameModes);
                
                spkPerFrm = cell(length(trialModes), nNrep); 
                for fi = 1:nNrep                    
                    
                    if any(strcmp(trialModes, 'all'))
                        spkPerFrm{1, fi} = mid_generateSpikeFile(Gid, cellId, windowMode_i, 'all', frameModes{fi}, responseType, responseToEachStim);
                    end
                    if any(strcmp(trialModes, 'odd'))
                        spkPerFrm{2, fi} = mid_generateSpikeFile(Gid, cellId, windowMode_i, 'odd', frameModes{fi}, responseType, responseToEachStim);
                    end
                    if any(strcmp(trialModes, 'even'))
                        spkPerFrm{3, fi} = mid_generateSpikeFile(Gid, cellId, windowMode_i, 'even', frameModes{fi}, responseType, responseToEachStim);
                    end
                    
                    if isequal(sort(trialModes), {'all', 'even', 'odd'})
                        diff = sum( (spkPerFrm{2, fi} + spkPerFrm{3, fi}) - spkPerFrm{1, fi} );
                        assert(abs(diff) < 1e-10);
                    end

                end
                
                if nNrep > 1
%                     for i = 1:nTrialModes
                    spf_1 = cell(1, nNrep);
                    
                    spf{1} = spf_1;
%                     spf{1} = nan;
                    for j = 2:nNrep
                        spf{j} = mid_generateSpikeFile(Gid, cellId, windowMode_i, frameModes{j}, 'all', relContrOfFrameToSpike);
                        diff = sum(spf{1})-sum(spf{j});
                        if ~isnan(diff)
                            assert(abs(diff) < 1e-10);
                        end
                    end
                end
                                
                
			end
			progressBar;
            3;
        end
        return;
    end


    % get spikes per frame info
    redo = 1;
    
    spkFileName = mid_getSpikeFileName(Gid, cellId, windowMode, trialMode, frameMode, responseType);    
    if exist(spkFileName, 'file') && ~redo
        spksPerFrame_exact = nan; 
        return;
    end

    if length(varargin) == 1
        relContrOfFrameToSpike = varargin{1};
        
    else        
        if length(varargin) == 2
            [spkWindow_ms, spkProfile] = varargin{:};
        else
            PSTHdata = getPSTHforCell(Gid, cellId);            
            spkWindow_ms  = PSTHdata.timeWindow_ms  ; %[30 60]; % defaultWindow
            spkProfile = PSTHdata.windowProfile;  % [1 1]; % defaultProfile
        end                                
%         relContrOfFrameToSpike = getParsedSpikes('frame',  Gid, cellId, spkWindow_ms, spkProfile );            
        relContrOfFrameToSpike = getSpkResponseToEachFrame(Gid, cellId, spkWindow_ms, spkProfile, 'stimOrder');
        
%         imagesc( mean(reshape(mean(relContrOfFrameToSpike,1), [36, 10, 4]),3) )
        
    end    
    stimType = getGratingStimType(Gid);        

%     if ~exist('windowMode', 'var')
%         windowMode = 'best';        
%     end
    
    if ~exist('frameMode', 'var')
        frameMode = 'all';        
    end
    if ~exist('trialMode', 'var')
        trialMode = 'all';
    end
    nTrials = stimType.nTrials;  
    
    
    if strcmp(frameMode, 'all') % considers all movies separately, even if identical
        nReps = nTrials;
    else
        nReps = sscanf(frameMode, '%drep');  
        assert(strcmp(frameMode, sprintf('%drep', nReps)));            
        assert(nReps <= nTrials);
    end
    if (nReps > 1) && ~strcmp(trialMode, 'all')
        error('Need nReps = 1 for odd/even trialMode');
    end
        
    switch trialMode
        case 'all', trial_idxs_use = 1:nTrials;
        case 'odd', trial_idxs_use = 1:2:nTrials;
        case 'even', trial_idxs_use = 2:2:nTrials;
    end
    nTrialsUse = length(trial_idxs_use);
    
    nTrialsPerRep = nTrialsUse/nReps;       assert(nTrialsPerRep == round(nTrialsPerRep));
    trial_idxs = reshape(trial_idxs_use, [nReps, nTrialsPerRep]);
    trial_idxs = mat2cell(trial_idxs, ones(nReps, 1), nTrialsPerRep);

    if iscell(relContrOfFrameToSpike)
        relContrOfFrameToSpike = [relContrOfFrameToSpike{:}];
    end
    frameStimIds = getStimulusFrameSequence(Gid);
    uStimIds = uniqueInOrder(frameStimIds); % keep original (pseudorandom) order    

%     stim_idx = ord(uStimIds);
    stim_idx = (uStimIds);

    relContrOfFrameToSpike = relContrOfFrameToSpike(:,stim_idx); % re-order spike response values to match (pseudorandom) order of stimulus
    
%     spksForStims = cell(1, nReps);
%     for rep_i = 1:nReps
%         rep_trial_idx = trial_idxs{rep_i}; % for the i_th repeat, are using these trial indices.
%         spksForStims{rep_i} = arrayfun(@(stim_i) sum( relContrOfFrameToSpike ( stimIdxs{stim_i}(rep_trial_idx) ) ), 1:nStims);
%     end
    
    spksForStims = cell(1, nReps);
    for rep_i = 1:nReps
        rep_trial_idxs = trial_idxs{rep_i}; % for the i_th repeat, are using these trial indices.        
        spksForStims{rep_i} = sum( relContrOfFrameToSpike(rep_trial_idxs, :), 1);        
    end            
    % reorder according to stimulus order.



    spksPerFrame_exact = [spksForStims{:}];  % concatenate
%     if strcmp(frameMode, 'all')
%         spksPerFrame_exact = relContrOfFrameToSpike;
%     end        
    
    uint8max = double(intmax('uint8'));
    spksPerFrame_uint8 = uint8(  spksPerFrame_exact/max(spksPerFrame_exact) * uint8max ); % rescale, convert to uint8
    
    doCheckIfFileExists = 1;
    if doCheckIfFileExists && exist(spkFileName, 'file')
        %%
        fileId_saved = fopen(spkFileName, 'r');
        spksPerFrame_uint8_saved = fread(fileId_saved, length(spksPerFrame_uint8),  'uint8');
        mtch = isequal(spksPerFrame_uint8_saved(:), spksPerFrame_uint8(:) );
        cc = corr( double(spksPerFrame_uint8_saved(:)), double(spksPerFrame_uint8(:)));
%         fprintf('%47s : match = %d. cc = %.3f\n', basename( spkFileName), mtch, cc);
        fclose(fileId_saved);
    else
    
        % open file    
        fileId = fopen(spkFileName, 'w');

        % save data to file
        fwrite(fileId, spksPerFrame_uint8, 'uint8');        
        fclose(fileId);  
        fprintf('Saved %47s\n', spkFileName);
    end    

end


%{
    spksPerFrame = cell(1, length(uMovieIds));
    for mi = 1:length(uMovieIds)
        m_idxs = movieIdxs{mi};
        spksPerFrame{mi} = relContrOfFrameToSpike{m_idxs(1)};
        for movieRep_i = 2:length(m_idxs)
            spksPerFrame{mi} = spksPerFrame{mi} + relContrOfFrameToSpike{m_idxs(movieRep_i)};
        end                    
    end

%}

%{
case 'uMovie'  % groups identical movies together
    movieIds = sd.stimulusInfo.movieIds;

    if length(movieIds) == 1
        spksPerFrame = [relContrOfFrameToSpike{:}];
    else
        [uMovieIds, movieIdxs] = uniqueListInOrder(movieIds);
        spksPerFrame = cellfun(@(idxs) sum( vertcat(relContrOfFrameToSpike{idxs}), 1), movieIdxs, 'un', 0);
        spksPerFrame = [spksPerFrame{:}];
    end
%}

%{    

%}

%     spksForStims = cell(1, nReps);
%     for rep_i = 1:nReps
%         rep_trial_idx = trial_idxs{rep_i}; % for the i_th repeat, are using these trial indices.
%         for stim_i = 1:nStims
%             spksForStims{rep_i}(stim_i) = sum( relContrOfFrameToSpike ( stimIdxs{stim_i}(rep_trial_idx) ) );
%         end
%     end          