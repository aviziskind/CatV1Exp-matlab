function varargout = getNoiseStimulusFrame(varargin)
    % with @getFrame = @getNoiseStimulusFrame;
    % first time: call with 
    %       getFrame('load', 'Gid', Gid);
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
    
    
    persistent hnd Did fileId;
    persistent nrows ncols pixelsPerFrame lastFrameRead nFramesEachPres cumulativeFrames;
    persistent pixelRange maxIntensity minIntensity Amax Amin;

    INT = 101;
    if (nargin == 0)

        frameToRead = lastFrameRead+1;

    
    elseif (nargin > 0)

        if isnumeric(varargin{1})
        
            frameToRead = varargin{1};            
        
        elseif ischar(varargin{1})
            switch varargin{1}
                case 'load' 
                    % Can either open with Gid/Did (frameId counts from beginning of all noise presentations shown
                    % to the group) or with noisePresId (counts from beginning of specific presentation)
                    hnd = dbOpenExpDb;
                    scaleMethod = 'normal';
                    tableName = 'TBL_NOISE_PRES';

                    switch varargin{2}
                        case {'Gid', 'Did'}
                            Did = dbLookup('Did',  varargin{2}, varargin{3});
                            criteria = {'DATAFILE_ID', Did};

                        case 'presId'
                            noisePresId = varargin{3};
                            criteria = {'NOISE_PRES_ID', noisePresId};                            
                    end
                    
                    [nSustainedFrames, nrows, ncols, nGradations, maxIntensity, minIntensity] = getFieldsFromDatabaseTable(hnd, ...
                        {'LNG_N_SUSTAINED_FRM', 'LNG_N_ROWS', 'LNG_N_COLUMNS', 'LNG_N_GRADATIONS', 'DBL_MAX_INTENSITY', 'DBL_MIN_INTENSITY'}, tableName, criteria);

                    Amax = maxIntensity;
                    Amin = minIntensity;
                    nrows = unique(nrows);
                    ncols = unique(ncols);
                    pixelsPerFrame = nrows * ncols;
                    nFramesEachPres = nSustainedFrames;
                    cumulativeFrames = [cumsum(nFramesEachPres)];
                    pixelRange = [0 unique(nGradations)-1];

                    noiseFilename = getName('noiseStimFile', 'Did', Did);
                    fileId = fopen(noiseFilename, 'r');                                        
                    
                    if nargin > 3
                        if strcmp(varargin{4}, 'scaling')
                            scaleMethod = varargin{5};
                        end
                    end
                    
                    if strcmp(scaleMethod, 'normal')
                        Amax = maxIntensity;
                        Amin = minIntensity;
                    elseif strcmp(scaleMethod, 'aroundZero')
                        d = mean([minIntensity, maxIntensity], 2);
                        Amax = maxIntensity - d;
                        Amin = minIntensity - d;
                    elseif strcmp(scaleMethod, 'int')
                        Amax = INT;
                    end                    
                    
                case 'size'
                    varargout = {[nrows ncols]};

                case 'nTotalPres'
                    varargout = {sum(nFramesEachPres)};

                case 'nFramesEachPres'
                    varargout = {nFramesEachPres};
                    
                case 'close'
                    fclose(fileId);
                    
            end
            return;
            
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% include pre/post blank frames????  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nidIndex = find(frameToRead <= cumulativeFrames, 1, 'first');
    if isempty(nidIndex)
        error('Invalid Noise Frame (frame Index Out of bounds)');
    end

    frameIdThisNoisePres = frameToRead - el([0; cumulativeFrames], nidIndex);
    lastFrameRead = frameIdThisNoisePres;
    
        fseek(fileId, (frameIdThisNoisePres-1)*pixelsPerFrame, 'bof');
        if (status == -1)
            error('stim:read', ferror(fileId);
        end

        frame = fread(fileId, [nrows, ncols] )';  % transpose because fread reads in column order, whereas the stimulus was presented in row order.
        if ~(Amax == INT)
            frame = Amin(midIndex) + (diff([Amin(midIndex) Amax(midIndex)])/diff(pixelRange)) * double(frame);
        end

    varargout = {frame};
    

end