function varargout = getMovieStimulusFrame(varargin)
    % with @getFrame = @getMovieStimulusFrame;
    % first time: call with 
    %       getFrame('load', Gid);
    % To retrieve frames with intensity adjusted so that the mean is 0, initialize with
    %       getFrame(..., 'scaling', 'aroundZero');
    % then just call with
    %       frame = getFrame(frameId);
    %  or:
    %       frame = getFrame;   % for the next frame after the one just retrieved
    % also available is:
    %       [nrows, ncols] = getFrame('size');
    % when done, please call
    %       getFrame('close');    
   
    persistent movieIds fileIds;
    persistent nrows ncols pixelsPerFrame nFramesEachPres cumulativeFrames0 lastFrameRead;
    persistent rescalePixels maxIntensity minIntensity Amax Amin;

    
    pixelRange = [0 255];     
    
    if (nargin == 0)

        frameToRead = lastFrameRead+1;

    
    elseif (nargin > 0)

        if isnumeric(varargin{1})
        
            frameToRead = varargin{1};            
        
        elseif ischar(varargin{1})
            switch varargin{1}
                case 'load'
                    % Can either open with Gid / Did (frameId counts from beginning of all movie presentations shown
                    % to the group) or with moviePresId (counts from beginning of specific presentation)
%                     hnd = dbOpenExpDb;
                    Gid = varargin{2};
                    
                    sd = siteDataFor(Gid);
                    stimInfo = sd.stimulusInfo;
                    
                    scaleMethod = 'normal';

                    [movieIds, uMovieFileNames, nFramesEachPres, nrows, ncols, maxIntensity, minIntensity] = deal(...
                        stimInfo.movieIds, stimInfo.movieFiles, stimInfo.nFramesPerPres, stimInfo.nrows, stimInfo.ncols, stimInfo.int_max, stimInfo.int_min);
                    
                    pixelsPerFrame = nrows * ncols;
                    nFramesEachPres_all = nFramesEachPres*ones(length(movieIds), 1);
                    cumulativeFrames0 = [0; cumsum(nFramesEachPres_all)];                    
                    
                    nPres = length(movieIds);
                    [uMovieIds, movieIdx] = uniqueList(movieIds);
                    
                    fgMoviePath = getName('fgMoviePath');
                    fileIds = zeros(1,nPres);                    
                    for mi = 1:length(uMovieIds)                        
                        movieFileName = [fgMoviePath, lower(uMovieFileNames{mi})];
                        [fileId_i, message] = fopen(movieFileName, 'r');  error(message);
                        fileIds(movieIdx{mi}) = fileId_i; 
                    end

                    if nargin > 3
                        if strcmp(varargin{3}, 'scaling')
                            scaleMethod = varargin{4};
                        end
                    end
                    
                    if strcmp(scaleMethod, 'normal')
                        rescalePixels = 1;
                        Amax = maxIntensity;
                        Amin = minIntensity;
                    elseif strcmp(scaleMethod, 'aroundZero') % image intensity from -1 to 1
                        rescalePixels = 1;
                        Amax = maxIntensity; 
                        Amin = -maxIntensity;
                        
%                         d = (maxIntensity-minIntensity);
%                         Amin = minIntensity - d;
%                         d = (minIntensity+maxIntensity)/2;
%                         Amax = maxIntensity - d;
%                         Amin = minIntensity - d;
                    elseif strcmp(scaleMethod, 'none')
                        rescalePixels = 0;
                    end

                case 'size'
                    if nargout == 1                        
                        varargout = {[nrows ncols]};
                    elseif nargout == 2
                        varargout = {nrows, ncols};
                    end

                case 'nTotalFrames'
                    varargout = {sum(nFramesEachPres)};

                case 'nFramesEachPres'
                    varargout = {nFramesEachPres};
                    
                case 'close'
                    [uMovieIds, firstMovieIdx] = unique(movieIds, 'first');
                    for m = 1:length(uMovieIds)
                        fclose(fileIds(firstMovieIdx(m)));
                    end

            end
            return;
            
        end
    end
        
    midIndex = find(frameToRead <= cumulativeFrames0, 1, 'first')-1;
    if isempty(midIndex)
        error('Invalid Movie Frame (frame Index Out of bounds)');
    end

    lastFrameRead = frameToRead;
    frameIdThisMovie = frameToRead - cumulativeFrames0(midIndex);
        
        fid = fileIds(midIndex);
        fseek(fid, (frameIdThisMovie-1)*pixelsPerFrame, 'bof');   error(ferror(fid));
        frame = fread(fid, [nrows, ncols])';  % transpose because fread reads in column order, whereas the stimulus was presented in row order.
        
        if rescalePixels            
            frame = Amin + (Amax-Amin)/diff(pixelRange) * double(frame);
        end
        
    varargout = {frame};
    

end


%     global checkForMissingFrame;

%     frameSequence = 1:cumulativeFrames0(end);                    
%     if checkForMissingFrame
%         [anyMissedFrames, actualFrameSequence] = getMissedFrames('Did', Did);    
%         if anyMissedFrames
%             frameSequence = actualFrameSequence;
%         end 
%     end

%     frameToReadIdx = frameSequence(frameToRead);
%     midIndex = find(frameToReadIdx <= cumulativeFrames0, 1, 'first')-1;
%     if isempty(midIndex)
%         error('Invalid Movie Frame (frame Index Out of bounds)');
%     end
% 
%     lastFrameRead = frameToRead;
%     frameIdThisMovie = frameToReadIdx - cumulativeFrames0(midIndex);



% 
%                     [movieIds, nrows, ncols, maxIntensity, minIntensity] = getFieldsFromDatabaseTable(hnd, ...
%                         {'MOVIE_ID', 'LNG_N_ROWS', 'LNG_N_COLUMNS', 'DBL_MAX_INTENSITY', 'DBL_MIN_INTENSITY'}, tableName, criteria);
% 
%                     [nSustainedFrames, nDisplayedFrames, preBlankFrames, postBlankFrames] = getFieldsFromDatabaseTable(hnd, ...
%                         {'LNG_N_SUSTAINED_FRM', 'LNG_N_DISPLAYED_FRM', 'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM'}, tableName, criteria);
% 
%                     if any ( nDisplayedFrames < nSustainedFrames + preBlankFrames + postBlankFrames)
%                         error('stim:interrupted', 'Movie was interrupted: Not all frames were displayed');
%                     end
