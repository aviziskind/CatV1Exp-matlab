function varargout = getMseqStimulusFrame(varargin)
    % with @getFrame = @getMseqStimulusFrame;
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
    
    
    persistent hnd Did fileId;
    persistent nrows ncols pixelsPerFrame lastFrameRead nFramesEachPres cumulativeFrames fileLength;
    persistent maxIntensity minIntensity Amax Amin tapRegister nMseqBits;

    INT = 101;
    if (nargin == 0)

        frameToRead = lastFrameRead+1;

    
    elseif (nargin > 0)

        if isnumeric(varargin{1})
        
            frameToRead = varargin{1};            
        
        elseif ischar(varargin{1})
            switch varargin{1}
                case 'load' 
                    % Can either open with Gid/Did (frameId counts from beginning of all mseq presentations shown
                    % to the group) or with mseqPresId (counts from beginning of specific presentation)
                    hnd = dbOpenExpDb;
                    scaleMethod = 'normal';
                    tableName = 'TBL_MSEQ_PRES';

                    switch varargin{2}
                        case {'Gid', 'Did'}
                            Did = dbLookup('Did',  varargin{2}, varargin{3});
                            criteria = {'DATAFILE_ID', Did};

                        case 'presId'
                            mseqPresId = varargin{3};
                            criteria = {'MSEQ_PRES_ID', mseqPresId};                            
                    end
                    
                    [tapRegister, nMseqBits, nDisplayedFrames, maxIntensity, minIntensity] = getFieldsFromDatabaseTable(hnd, ...
                        {'LNG_TAP_REGISTER', 'LNG_MSEQ_BITS', 'LNG_N_DISPLAYED_FRM', 'DBL_MAX_INTENSITY', 'DBL_MIN_INTENSITY'}, tableName, criteria);

                    nrows = 64; % ???? not sure what the screen dimensions were...
                    ncols = 64; 
                    Amax = maxIntensity;
                    Amin = minIntensity;
                    pixelsPerFrame = nrows * ncols;
                    nFramesEachPres = nDisplayedFrames;
                    cumulativeFrames = [cumsum(nFramesEachPres)];

                    mseqFilename = getName('mseqStimFile', tapRegister, nMseqBits );
                    fileId = fopen(mseqFilename, 'r');                                        
                    
                    fseek(fileId, 0, 'eof');
                    fileLength = ftell(fileId);
                    
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

    midIndex = find(frameToRead <= cumulativeFrames, 1, 'first');
    if isempty(midIndex)
        error('Invalid Mseq Frame (frame Index Out of bounds)');
    end

    frameIdThisMseqPres = frameToRead - el([0; cumulativeFrames], midIndex);
    lastFrameRead = frameIdThisMseqPres;

        % Read frame from file (handle looping to beginning of file when reach end of file)
        posStart = mod( (frameIdThisMseqPres-1)*pixelsPerFrame, fileLength);
        fseek(fileId, posStart, 'bof');
        
        if (fileLength-posStart) >= pixelsPerFrame
            frameData = fread(fileId, pixelsPerFrame);
        else
            L = [(fileLength-posStart), pixelsPerFrame-(fileLength-posStart)]; % loop back to beginning
            frameData = fread(fileId, L(1)); 
            fseek(fileId, 0, 'bof');
            frameData = [frameData; fread(fileId, L(2))];
        end
            
        frame = reshape(frameData, [nrows, ncols] )';  % transpose because fread reads in column order, whereas the stimulus was presented in row order.
        if ~(Amax == INT)
            frame = Amin(midIndex) + (diff([Amin(midIndex) Amax(midIndex)])) * double(frame);
        end

    varargout = {frame};
    

end