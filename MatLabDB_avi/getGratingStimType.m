function varargout = getGratingStimType(arg)

%     persistent allFgGids allFgStimTypes
    if isnumeric(arg)
        Gid = arg;
        siteData = siteDataFor(Gid);
        stimType = siteData.stimType;
    elseif ischar(arg)
        stimType = arg;
    end
        
    % usually don't want the first 2 classifications.
    if isnan( str2double(stimType(1)) )  % ie. first character is a letter (not a number)
%         idx_colon = strfind(stimType, ':');
        [str1, remain1] = strtok(stimType, ':');
        [str2, remain2] = strtok(remain1(2:end), ':');
        
        stimType = remain2(2:end);
    end

    A = sscanf(stimType, '%dx%dx%d(%dx%d)');
    [nOri, nSp, nPh, nUniqueSeq, nSeqRep] = dealV(A);
    
    gratingType = iff(strcmp(str2, 'Flashed_Gratings'), 'flashed', 'drifting');
    isCphFlashed = ~isempty(strfind(stimType, 'C'));
    if strcmp(gratingType, 'drifting')
        driftingSubtype = str2;
    else
        driftingSubtype = '';
    end
        
    if nargout <= 1
        s = struct('nOri', nOri, 'nSpf', nSp, 'nPh', nPh, 'nTrials', nUniqueSeq*nSeqRep, 'nUniqueSeq', nUniqueSeq, 'nSeqRep', nSeqRep, ...
                'gratingType', gratingType, 'isCphFlashed', isCphFlashed, 'driftingType', driftingSubtype);
        
        varargout = {s};                
    else
    
        varargout = {nOri, nSp, nPh, nUniqueSeq*nSeqRep, isCphFlashed};
    end       
    
end



%     if isempty(allFgGids)
%         S_mcg = load('cellsGroups_movie_fg');
%         mcg = S_mcg.movieGroups_fg;
%         allFgGids = [mcg.Gid];
%         allFgFrmLengths = {mcg.frameLength_ms};
%         allFgStimTypes = {mcg.stimType};
%         n = length('Movie:Flashed_Gratings:');        
%         allFgStimTypes = cellfun(@(s, f) [s(n+1:end) '_' num2str(round(f))], allFgStimTypes, allFgFrmLengths, 'un', 0);        
%     end
    
%     Gid_idx = find(allFgGids == Gid, 1);
