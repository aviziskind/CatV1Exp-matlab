function spikeFileName = mid_getSpikeFileName(Gid, cellId, timeWindow, trialMode, frameMode, responseType)
    sd = siteDataFor(Gid);
    allMovieIds = sd.stimulusInfo.movieIds;    
    nStim = length(sd.ori_deg)*length(sd.spPeriod_pix)*length(sd.spPh_deg);
    nRepsPerMovie = sd.stimulusInfo.nFramesPerPres / nStim;
    nRepsTotal = length(allMovieIds) * nRepsPerMovie;

    if ~exist('timeWindow', 'var') || isempty(timeWindow) || strcmp(timeWindow, 'best') % 'best', or specific window (30-60, 60-90)
        timeWindow_str = '';
    elseif isnumeric(timeWindow) && length(timeWindow) == 2        
        timeWindow_str = sprintf('__%d-%d', floor(timeWindow(1)), floor(timeWindow(2)));
    end
        
    if ~exist('trialMode', 'var') || isempty(trialMode) || strcmp(trialMode, 'all')  % 'all', 'odd', or 'even'
        trialMode_str = '';
    else 
        trialMode_str = ['__' trialMode];
    end
        
    if ~exist('frameMode', 'var') || isempty(frameMode)
        frameMode = '1rep';
    end
    
    if ~exist('responseType', 'var') || isempty(responseType) || strcmp(responseType, 'raw')
        responseType_str = '';
    elseif strcmp(responseType, 'gainCorrected')
        responseType_str = '_GC';
    end
    
    if strcmp(frameMode, 'all')
        nReps = nRepsTotal;
    else
        nReps = sscanf(frameMode, '%drep');
        assert(strcmp(frameMode, sprintf('%drep', nReps)));
        assert(nReps <= nRepsTotal);
        nTrialsPerRep = nRepsTotal/nReps;
        if nTrialsPerRep ~= round(nTrialsPerRep)
            error('nReps must be a divisor of the total number of Trials');
        end        
    end
    
    
    file_path = getName('compiledSpikesPath');        
    filename = sprintf('Group_%d__cell_%d__%d_rep%s%s%s', Gid, cellId, nReps, timeWindow_str, trialMode_str, responseType_str);    
    
    spikeFileName = [file_path filename '.spk'];

end

