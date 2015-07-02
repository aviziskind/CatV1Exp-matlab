function frameDims = mid_generateStimulusFile(Gid, frameMode, movieFileName)

    % Open original stimulus file
    redo = 0;    
    if ~exist('frameMode', 'var')
        frameMode = 'all';
    end    
    if nargin < 3
        movieFileName = mid_getMovieStimFileName(Gid, frameMode);
    end
    if exist(movieFileName, 'file') && ~redo
        return;
    end
    
    sd = siteDataFor(Gid);    
    stimType = getGratingStimType(sd.stimType);
    nRepsTotal = stimType.nTrials;    
    
    if strcmp(frameMode, 'all') 
        nReps = nRepsTotal;
    else
        nReps = sscanf(frameMode, '%drep');  
        assert(strcmp(frameMode, sprintf('%drep', nReps)));            
        assert(nReps <= nRepsTotal);
    end
    
    frameIdSequence = getStimulusFrameSequence(Gid, 'OSP');
    uStimIdsInOrder = uniqueInOrder(frameIdSequence);
    [uStimIds, idxFirstOccurence] = unique(frameIdSequence, 'first');
%     stim_idx = ord(uStimIds);
    
    frameOrder = idxFirstOccurence(uStimIdsInOrder);
    frameIds = repmat(frameOrder(:)', [1, nReps]);

%     if strcmp(frameMode, 'all') 
%         frameIds = 1:length(frameIdSequence);
%     end
    
    % Open original stimulus files
    getFrame = getFrameRetrieverFunction(Gid);
    getFrame('load', Gid, 'scaling', 'none');
    frameDims = getFrame('size');
    
    % Open new stimulus file;    
    fileId = fopen(movieFileName, 'w');
    fprintf('Creating %s\n', movieFileName);
    % Copy from original --> new stimulus file
    for fi = 1:length(frameIds)
        frame = getFrame(frameIds(fi));
        fwrite(fileId, frame, 'uint8');        
    end
    
    % close stimulus files
    fclose(fileId);  
    getFrame('close');

end

%{
    switch frameMode
        case 'uStim'  % most compact: groups the stimuli according to type (only for Flashed Grating stimuli)           
            [frameIdSequence] = getStimulusFrameSequence(Gid, 'OSP');
            [uStimIds, idx_firstOccurence] = unique(frameIdSequence, 'first');
            frameIds = idx_firstOccurence;
            
        case 'uMovie'  % groups identical movies together        
            all_nFramesEachPres = nFramesEachPres(ones(length(movieIds), 1));
            
            uMovieIds = uniqueListInOrder(movieIds);
            frameIds = cell(1, length(uMovieIds));
                
            movieFrameIds = [0; cumsum(all_nFramesEachPres)];
            for mi = 1:length(uMovieIds)
                movieInd = find(movieIds == uMovieIds(mi), 1 );                    
                frameIds{mi} = movieFrameIds(movieInd)+1:movieFrameIds(movieInd+1) ; 
            end
            frameIds = [frameIds{:}];
            
        case 'all' % considers all movies separately, even if identical
            nTotalFrames = nFramesEachPres * length(movieIds);            
            frameIds = 1:nTotalFrames;
    end
%}


%{    
    gids = getAllGids('f'); n = length(gids);
    for i = 1:length(gids)
        mid_generateStimulusFile(gids(i), 'all');
        mid_generateStimulusFile(gids(i), '2rep');
        mid_generateStimulusFile(gids(i), '1rep');
    end
%}