function varargout = getGratingStimulusFrame(varargin)
    % with @getFrame = @getGratingStimulusFrame;
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
    
    persistent hnd Did;
    persistent lastFrameMade;
    persistent isStaticGrating isSquareGrating cumulativeFrames;
    persistent fps ori_deg spatPeriod_pix spatPhase_deg tempPeriod_sec tempPhase_deg; %theta omega phi k 
    persistent nFramesEachPres pixelRange maxIntensity minIntensity Amax Amin;

    INT = 101;
    if (nargin == 0)

        frameToMake = lastFrameMade+1;

    
    elseif (nargin > 0)

        if isnumeric(varargin{1})
        
            frameToMake = varargin{1};
        
        elseif ischar(varargin{1})
            switch varargin{1}
                case 'load' 
                    % Can either open with Gid (frameId counts from beginning of all noise presentations shown
                    % to the group) or with gratingPresId (counts from beginning of specific presentation)
                    hnd = dbOpenExpDb;
                    tableName = 'TBL_GRATING_PRES';
                    scaleMethod = 'normal';
                    
                    switch varargin{2}
                        case {'Gid', 'Did'}
                            Did = dbLookup('Did',  varargin{2}, varargin{3});
                            criteria = {'DATAFILE_ID', Did};

                        case 'presId'
                            gratingPresId = varargin{3};
                            criteria = {'GRATING_PRES_ID', gratingPresId};
                            Did = unique( getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', tableName, criteria) );
                    end
                   
                    fieldsToGet = {'DBL_ORIENTATION_DEGR', 'DBL_SPATIAL_PERIOD_PIX', 'DBL_SPATIAL_PHASE_DEGR', ...
                        'DBL_TEMP_PERIOD_FRM',  'DBL_TEMP_PHASE_DEG', ...
                        'BLN_STATIC_GRATING', 'BLN_SQUARE_GRATING', 'LNG_N_SUSTAINED_FRM', 'DBL_MAX_INTENSITY', 'DBL_MIN_INTENSITY' };

                    [ori_deg, spatPeriod_pix, spatPhase_deg, tempPeriod_frm, tempPhase_deg, isStaticGrating, isSquareGrating, nFramesEachPres, maxIntensity, minIntensity] = ...
                        getFieldsFromDatabaseTable(hnd, fieldsToGet, tableName, criteria);

                    isStaticGrating = (isStaticGrating == -1) || (tempPeriod_frm >= 50000);  % ie. 0 is false, -1 is true.
                    isSquareGrating = (isSquareGrating == -1);  % ie. 0 is false, -1 is true.
%                     k = 2*pi ./ spatPeriod_pix;
%                     phi = deg2rad(spatPhase_deg);
%                     theta = deg2rad(ori_deg);
                    if isStaticGrating
                        tempPeriod_sec = inf;
%                         omega = (2*pi) ./ period_frm;
                    else
                        tempPeriod_sec = tempPeriod_frm / fps;
%                         omega = 0;
                    end
                    cumulativeFrames = [cumsum(nFramesEachPres)];
                    fps = dbConvertTimeMeasures(Did, 1, 'sec', 'frame');
                    varargout = {sum(nFramesEachPres)};

                    pixelRange = [-1 1];
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
                    varargout = {size(getGratingStimulusFrame(1)) };

                case 'nTotalPres'
                    varargout = {sum(nFramesEachPres)};

                case 'nFramesEachPres'
                    varargout = {nFramesEachPres};
                    
                    
                case 'close'

                    
            end
            return;
            
        end
    end

    gi = find(frameToMake <= cumulativeFrames, 1, 'first');
    if isempty(gi)
        error('Invalid Grating Frame (frame Index Out of bounds)');
    end

    lastFrameMade = frameToMake;
    frameIdThisGratingPres = frameToMake - el([0; cumulativeFrames], gi);
    
%                     k = 2*pi ./ spatPeriod_pix;
%                     phi = deg2rad(spatPhase_deg);
%                     theta = deg2rad(ori_deg);
    
    frameDims = [1024, 768]; % full field grating

        t_sec = frameIdThisGratingPres / fps;
        frame = generateGratingFrame(frameDims, ori_deg(gi), spatPeriod_pix(gi), spatPhase_deg(gi), tempPeriod_sec(gi), tempPhase_deg, t_sec, isSquareGrating(gi) );
                    
        frame = Amin(gi) + (diff([Amin(gi) Amax(gi)])/diff(pixelRange)) * frame;
        if ~(Amax == INT)
            frame = Amin(midIndex) + (diff([Amin(midIndex) Amax(midIndex)])/diff(pixelRange)) * double(frame);
        end
        
    varargout = {frame};
            
end


