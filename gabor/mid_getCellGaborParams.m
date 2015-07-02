function [gparams, rsqr, MID_fit, s] = mid_getCellGaborParams(Gid, cellId, timeWindow, trialMode, responseType, jackknifeIdx, data_input, calcMID_fit_arg)

    persistent allCellGaborParams saveCount 
    
    cellGaborParams_file = [CatV1Path 'MatLabDB_avi' filesep 'allCellGaborParams.mat'];
        
    redo_all = 0;
    redo_current = 0;
    saveCountSpacing = 1;
    redo_ifLowerDsampAvailable = 1;
    haveInputs = nargin > 5 && ~isempty(data_input);
    calcMID_fit = exist('calcMID_fit_arg', 'var') && isequal(calcMID_fit_arg, 1);
    returnEmptyIfDontHaveData = 0;
    
    if strcmp(Gid, 'fix')
        if isempty(allCellGaborParams)
            S_file = load(cellGaborParams_file);
            allCellGaborParams = S_file.allCellGaborParams;
        end
        fld_names = fieldnames(allCellGaborParams);
        for i = 1:length(fld_names)
            [Gid, cellId, timeWindow, trialMode, responseType, jackIdx] = parseFieldName(fld_names{i}, 1);
            [gparams, rsqr] = mid_getCellGaborParams(Gid, cellId, timeWindow, trialMode, responseType, jackIdx);            
        end        
        saveCount = 1;
        mid_getCellGaborParams('save');        
    end
    
    if strcmp(Gid, 'save')
        %%
        if ~isempty(saveCount) && saveCount > 0
            %%
            fprintf('[saved %d gabor parameter sets]\n', length(fieldnames(allCellGaborParams)))
            allCellGaborParams = orderfields(allCellGaborParams);
            save(cellGaborParams_file, 'allCellGaborParams', '-v6');                
            saveCount = 0;
        end
        %%
        return;
    elseif any(strcmp(Gid, {'used', 'all'}))
        if strcmp(Gid, 'used')
            S = load('usedCells.mat');
            Gids = S.usedGids;
            cellIds = S.usedCellIds;            
        elseif strcmp(Gid, 'all')
            S = load('allCells.mat');
            Gids = S.allGids;
            cellIds = S.allCellIds;            
        end                        
        for i = 1:length(cellIds)
            Gid = Gids(i); cellId = cellIds(i);
            fprintf('(%d) Gid = %d, cellId = %d : ', i, Gid, cellId);
            [gparams, rsqr, MID_fit, s] = mid_getCellGaborParams(Gid, cellId);
            fprintf('Done : (dsamp = %d) rsqr = %.2f\n\n', s.downSampFactor, rsqr);
            3;
        end
        mid_getCellGaborParams('save');
        return;
    end

    
    if isempty(allCellGaborParams)        
        if exist(cellGaborParams_file, 'file') && ~redo_all
            S_file = load(cellGaborParams_file);
            allCellGaborParams = S_file.allCellGaborParams;
        else
            allCellGaborParams = struct;
        end        
        saveCount = 0;        
    end
    
    cell_fld_name = getGaborParamsFieldName(Gid, cellId, timeWindow, trialMode, responseType, jackknifeIdx, 1);    
    
    if ~haveInputs    
        redo_now = 0;
        if isfield(allCellGaborParams, cell_fld_name)
            s = allCellGaborParams.(cell_fld_name);
            if isfield(s, 'MID')
                MID_fit = getMIDfitToGabor(Gid, s.params);
                assert(isequal(s.MID_fit, MID_fit));
                
                s = rmfield(s, 'MID');
                s = rmfield(s, 'MID_fit');                
                allCellGaborParams.(cell_fld_name) = s;
            end            
            if isfield(s, 'all_p_est')
                s = rmfield(s, 'all_p_est');
                allCellGaborParams.(cell_fld_name) = s;
            end
            
            mid_fileName = mid_getPreferredMIDfile(Gid, cellId, timeWindow, trialMode, responseType);        
            if redo_ifLowerDsampAvailable && (s.downSampFactor > 1) && ~isempty(strfind(mid_fileName, '_d1.mat'))
                redo_now = 1;
            end
        end

        if (~isfield(allCellGaborParams, cell_fld_name) || redo_current || redo_now)
            if returnEmptyIfDontHaveData
                [gparams, rsqr, MID_fit, s] = deal([]);
                return;
            else
                s = calcGaborParameters(Gid, cellId, timeWindow, trialMode, responseType, jackknifeIdx);        

                allCellGaborParams.(cell_fld_name) = s;
                saveCount = saveCount + 1;

                if saveCount >= saveCountSpacing
                    allCellGaborParams = orderfields(allCellGaborParams);
                    save(cellGaborParams_file, 'allCellGaborParams', '-v6');        
                    fprintf('[saved %d gabor parameter sets]\n', length(fieldnames(allCellGaborParams)))
                    saveCount = 0;
                end                
            end
        end
        
    elseif haveInputs
        gparams_input = data_input.gparams;
        MIDfit_input = data_input.MIDfit;
        rsqr_input = data_input.rsqr;
        
        s = getGaborParamStruct(Gid, cellId, timeWindow, trialMode, responseType, gparams_input, MIDfit_input, rsqr_input);       
        allCellGaborParams.(cell_fld_name) = s;
        mid_getCellGaborParams('save');
    end
        
        
    if nargout > 0
        s = allCellGaborParams.(cell_fld_name);
        gparams = s.params;
        rsqr = s.rsqr;
        MID_fit = [];
        if nargout >= 3 && calcMID_fit
            MID_fit = getMIDfitToGabor(Gid, gparams);
%             MID_fit = s.MID_fit;           
        end
            
        
    end


end


function s = getGaborParamStruct(Gid, cellId, timeWindow, trialMode, responseType, gparams, MID_fit, rsqr)

    [~, downSampFactor, frameMode] = mid_getPreferredMIDfile(Gid, cellId, timeWindow, trialMode, responseType);
%     [mid_fileName, downSampFactor, frameMode] = mid_getPreferredMIDfile(Gid, cellId, timeWindow, trialMode);
    
%     S = load(mid_fileName);
%     MID = S.MID;
    
    sd = siteDataFor(Gid);
    [xs, ys] = getStimulusXY(sd.stimulusInfo, downSampFactor);    
    
%     gparams = calculateGaborParamsFromZ(xs, ys, MID);
        
%     [XX, YY] = meshgrid(xs, ys);
%     Z = reshape(gaborP(gparams, [XX(:), YY(:)]), size(XX));
%     [p_fits{est_i}, resid] = nlinfit(XX, YY, @gaborP, p0, statset('maxiter', 200));
    
    s.frameMode = frameMode;
    s.downSampFactor = downSampFactor; 
    s.xs = xs;
    s.ys = ys;
%     s.MID = MID;
    s.all_p_est = gparams;
    s.params = gparams;
%     s.MID_fit = MID_fit;
    s.rsqr = rsqr;


end

function s = calcGaborParameters(Gid, cellId, timeWindow, trialMode, responseType, jackknifeIdx)
    fprintf('Calculating gabor fit for Gid = %d, cellId = %d, timewindow = %s, trialMode = %s, jackKnife = %s\n', ...
        Gid, cellId, num2str(timeWindow), trialMode, iff(isempty(jackknifeIdx), '[ALL TRIALS]', num2str(jackknifeIdx)) );
    [mid_fileName, downSampFactor, frameMode] = mid_getPreferredMIDfile(Gid, cellId, timeWindow, trialMode, responseType);
    
    S = load(mid_fileName);
    if isempty(jackknifeIdx)
        MID = S.MID;
    else
        MID = S.v_MID{1}(:,:,jackknifeIdx);
    end
    
%     [nx, ny] = size(MID);
    
    sd = siteDataFor(Gid);
    [xs, ys] = getStimulusXY(sd.stimulusInfo, downSampFactor);    
    
%     gparams = calculateGaborParamsFromZ(xs, ys, MID);
    if strcmp(getenv('computername'), 'AVI-PC') 
        nMaxFits = 250;
    else
        nMaxFits = 10000;
    end
    nMaxFits = 500;
    
    
    smooth_Ws = [0 1 2 3];
%     nMaxFits = 30;
    circ_flag = 0;
    
    nSmooth = length(smooth_Ws);
    all_p0 = cell(1, nSmooth);
    for sm_i = 1:nSmooth
        smoothW = smooth_Ws(sm_i);
        MID_sm = gaussSmooth(gaussSmooth(MID, smoothW , 1, circ_flag), smoothW, 2, circ_flag);
        fprintf('smoothing = %.1f (%d/%d)', smoothW, sm_i, nSmooth);
        all_p0{sm_i} =  estimateGaborParameters(xs, ys, MID_sm, 1);
    end
    all_p0 = cat(1, all_p0{:});
    [p_fit, rsqr, MID_fit, ~, ~, all_p0_ordered] = fitGaborParams(xs, ys, MID, all_p0, nMaxFits);        
        
%     all_p0 = all_p0(best_fit_order, :);
  
    s.frameMode = frameMode;
    s.downSampFactor = downSampFactor; 
    s.xs = xs;
    s.ys = ys;
%     s.MID = MID;
    s.all_p_est = all_p0_ordered;
    s.params = p_fit;
%     s.MID_fit = MID_fit;
    s.rsqr = rsqr;
    s.fitSettings = struct('smooth_Ws', smooth_Ws, 'nMaxFits', nMaxFits);
    
end


function [Gid, cellId, timeWindow, trialMode, responseType, jackknifeIdx] = parseFieldName(fld_name, chk_flag)
% [Gid_test, cellId_test, timeWindow_test, trialMode_test, responseType_test, jackIdx_test] = parseFieldName(fld_name);
    if ~isempty(strfind(fld_name, '__even'))
        trialMode = 'even';                
    elseif ~isempty(strfind(fld_name, '__odd'))
        trialMode = 'odd';                
    else
        trialMode = 'all';
    end

    if ~isempty(strfind(fld_name, '__29_62'))
        timeWindow = [29 62];                
    elseif ~isempty(strfind(fld_name, '__58_91'))
        timeWindow = [58 91];                
    else
        timeWindow = 'best';
    end            

    if ~isempty(strfind(fld_name, '__GC'))
        responseType = 'gainCorrected';                
    else 
        responseType = 'raw';                
    end
    
    
    jackknifeIdx = [];
    i_jack = strfind(fld_name, '__jack1');
    if ~isempty(i_jack)
        jackknifeIdx = str2double( fld_name(i_jack+6:end) );
    end
    
    A = sscanf(fld_name, 'GaborParams_Gid_%d__cell%d');
    Gid = A(1); cellId = A(2);
    
    if nargin > 1 && isequal(chk_flag, 1)
        fld_name_chk = getGaborParamsFieldName(Gid, cellId, timeWindow, trialMode, responseType);
        assert(isequal(fld_name_chk, fld_name));
    end
    
end

function fld_name = getGaborParamsFieldName(Gid, cellId, timeWindow, trialMode, responseType, jackIdx, chk_flag)

    trialMode_str  = iff( strcmp(trialMode, 'all'),                          '', sprintf('__%s',    trialMode) );
    timeWindow_str = iff( strcmp(timeWindow, 'best') || isempty(timeWindow), '', sprintf('__%d_%d', timeWindow));
    jackIdx_str = iff( ~isempty(jackIdx), sprintf('_jack%d', jackIdx), '');
    responseType_str = iff( strcmp(responseType, 'gainCorrected'), '__GC', '') ;
    
    fld_name = sprintf('GaborParams_Gid_%d__cell%d%s%s%s%s', Gid, cellId, timeWindow_str, trialMode_str, responseType_str, jackIdx_str);

    
    if nargin >= 7 && isequal(chk_flag, 1)
        
        [Gid_test, cellId_test, timeWindow_test, trialMode_test, responseType_test, jackIdx_test] = parseFieldName(fld_name);
        assert(isequal(Gid, Gid_test));
        assert(isequal(cellId, cellId_test));
        assert(isequal(timeWindow, timeWindow_test));
        assert(isequal(trialMode, trialMode_test));        
        assert(isequal(responseType, responseType_test));
        assert(isequal(jackIdx, jackIdx_test));
    end

end


%{
    S = load('usedCells.mat');
    Gids = S.usedGids;
    cellIds = S.usedCellIds;
    for i = 1:length(cellIds)
        Gid = Gids(i); cellId = cellIds(i);
        [gparams, rsqr, MID_fit] = mid_getCellGaborParams(Gid, cellId);
        fprintf('(%d) Gid = %d, cellId = %d, rsqr = %.2f\n', i, Gid, cellId, rsqr)
    end
%}