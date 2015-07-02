function varargout = getOspDataForPsthWindow(Gid, cellId, allBins, psthVals, l_bin, r_bin, windowProfile, dataToCalculate, redo_flag)

    persistent allWindOSPs saveCount;
    global psthStatsSettings

    if nargin > 1
%         fprintf('Get OSP Data Called with Gid=%d, cellId=%d, L=%d, R=%d, length(wind) = %d, data=%s\n', ...
%                 Gid, cellId, l_bin, r_bin, length(windowProfile), tostring(dataToCalculate) )
    end
    
    
    if isempty(psthStatsSettings)
        setPsthGlobalSettings;
    end
    
    ignoreFile = false;        
    redo_all = false;
    redo_current_cell = 0;
    redo_current_cell_data = 0;
%     saveCountSpacing = 200;
    saveCountSpacing = 300;
    checkClustIdsMatchSaved = false;

    removeDeletedInfoFromAllFields = true;
    
    redo_old_cell_data = 0;
    
    if isempty(saveCount)
        saveCount = 0;
    end           
    if exist('redo_flag', 'var') && ~isempty(redo_flag)
        if redo_flag == 1  %% intepret as: 'always redo'
            redo_current_cell_data = 1;
        elseif redo_flag > 1  %% interpret as date.
            redo_old_cell_data = 1;
            redo_cell_data_if_before = redo_flag;
        end
    end
        
%     redo_current_cell_data_now = redo_current_cell_data || (exist('redo_flag', 'var') && ~isempty(redo_flag));

    %     redo_current_stat_if_short = false;        
    if ischar(Gid) && strcmp(Gid, 'save')
        if isempty(allWindOSPs) || (saveCount == 0)
            return;
        end        
        fn = fieldnames(allWindOSPs);
        isFlashed  = ~isempty(strfind(fn{1}, '_f_'));
        isDrifting = ~isempty(strfind(fn{1}, '_d_'));
        assert(xor(isFlashed, isDrifting));
        gratingType = iff(isFlashed, 'flashed', 'drifting');        
    elseif isnumeric(Gid)        
        gratingType = flashedOrDrifting(Gid, 'str');  
    end   
        
    methodName = 'stimulus';
    methodName = methodName(1:2);    
    
    ospFcn = psthStatsSettings.ospPhCompressFcn;

    responseType = curResponseType('');
    gainCorr_str = switchh(responseType, {'raw', 'gainCorrected'}, {'', '_GC'});
    
    timeWindow_str = curTimeWindow('');
    if strcmp(gratingType, 'drifting') || strcmp(timeWindow_str, 'best')
        timeWindow_str = '';
    end
    
    cellClustStr = curGroupingType('');    
    fet_str = iff(curMatchDB, '_DB', ['_' curSortingFeatures('')]);
    
    dataMethod = [cellClustStr fet_str '_' gratingType(1) '_' methodName '_' ospFcn gainCorr_str timeWindow_str];
            
    cellFieldName_fcn = @(Gid, cellId) sprintf('Gid_%04d_cell_%s_%s', Gid, iff(cellId == -1, 'n1', num2str(cellId)), methodName);
    ospWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'allWindowOspData_' dataMethod '.mat'];            
    saveFile = false;
    %%
    dataCanCalculate = {'osp', 'osp_oe',  'osp_ph', 'osp_ph_oe',  'osp_ph_hoe',  'osp_full', 'osp_full_oe',  'osp_pred',  'osp_ph_jackTrials', 'osp_ph_oe_jackTrials',  'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm',   ...
        'phaseTC_CCs_oe', 'phaseTC_CC_ps_oe', 'phaseTC_CCs_hoe', 'phaseTC_CC_ps_hoe', 'phaseTC_CCs_fs', 'phaseTC_CC_ps_fs', 'phaseTC_Dots', 'F1oDC', ...
        'stimF1oDCs', 'stimF1oDC_jackStds',  'STA', 'STA_oe', 'STA_jack', 'STA_odd_jack', 'STA_even_jack', 'meanRate', 'tuningStats', 'windowStats'};
    
    dataToDelete = {'osp_odd', 'osp_even', 'osp_odd_ph', 'osp_even_ph', 'STA_odd', 'STA_even', 'phaseTC_CCs', 'phaseTC_CC_ps', 'phaseTC_CC_hoe', 'STA_oe_jack', ...
        'osp_gainCorrected', 'osp_gainCorrected_oe', 'osp_gainCorrected_hoe', ...
        ...'osp_ph_oe_jackTrials',  'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm', 'STA_jack', 'STA_odd_jack', 'STA_even_jack'
        };
    
%     dataCanCalculate(strCcmp(dataCanCalculate, dataToDelete)) = [];

   %%
    
    
    if ischar(Gid) && strcmp(Gid, 'save') && ~isempty(allWindOSPs)
        %%
        fprintf('[saving OSP Data file '); tic;
        fn = fieldnames(allWindOSPs);
%         fn = fn(1); % if do more than just 1 - have to adjust save file name
        for i = 1:length(fn)
            allCellData = allWindOSPs.(fn{i}); %#ok<NASGU>
            fprintf('(%s) ...', fn{i});
            ospWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'allWindowOspData_' fn{i} '.mat'];            
            save(ospWindowDataFile, '-struct', 'allCellData', '-v6');                            
        end
        t_sec = toc; fprintf('done (%.2f s)]\n', t_sec);                
        saveCount = 0;
        return;
    end
    
        
%     if ~exist('dataToCalculate', 'var') || isempty(dataToCalculate)
%         dataToCalculate = {'osp', 'osp_odd', 'osp_even'};
%     end
    if ischar(dataToCalculate)
        dataToCalculate = {dataToCalculate};
    end
%     redo_current_cell_data_now = any(strCcmp({'STA'}, dataToCalculate));
    
    if ~all( strCcmp(dataToCalculate, dataCanCalculate))
        error('Provided unknown variables to calculate');
    end
    
    if isempty(allWindOSPs)
        allWindOSPs = struct;
    end
    
    if ~isfield(allWindOSPs, dataMethod) || redo_all || ignoreFile
        if ~exist(ospWindowDataFile, 'file') || redo_all || ignoreFile
            allWindOSPs.(dataMethod) = struct;
        else    
            fprintf('Loading saved OSP Data (%s) ', basename(ospWindowDataFile) ); tic; 
            allWindOSPs.(dataMethod) = load(ospWindowDataFile);
            fprintf(' done. '); toc;
        end
    end

    if isempty(saveCount)
        saveCount = 0;
    end

     
    if removeDeletedInfoFromAllFields && 0
        %%
        fn = fieldnames(allWindOSPs);
        
        for i = 1:length(fn)
            allCellData = allWindOSPs.(fn{i}); 
            cell_field_names = fieldnames(allCellData);
            for cell_j = 1:length(cell_field_names)
                fprintf('%d/%d : %s', cell_j, length(cell_field_names), cell_field_names{cell_j})
                cell_fields = fieldnames( allCellData.(cell_field_names{cell_j}) );
                if any(strCcmp(cell_fields, dataToDelete))
                    3;
                    fprintf( '[%s] ', cellstr2csslist( intersect(cell_fields, dataToDelete)));
                    allCellData.(cell_field_names{cell_j}) = rmfields( allCellData.(cell_field_names{cell_j}), ...
                        dataToDelete);
                end
                fprintf('\n')
            end
            allWindOSPs.(fn{i}) = allCellData;                
        end
        t_sec = toc; fprintf('done (%.2f s)]\n', t_sec);                
        
        
    end
        
    % if doesn't exist in loaded file, calculate it
    cellFieldName = cellFieldName_fcn(Gid, cellId);    
    haveSomeDataForThisCell = isfield(allWindOSPs.(dataMethod), cellFieldName);

    if ~curMatchDB && strcmp(cellClustStr, 'cells') && checkClustIdsMatchSaved
        curClustIds = getClusterCellAssignment(Gid, cellId);
        
        if haveSomeDataForThisCell && ~curMatchDB && strcmp(cellClustStr, 'cells')
            savedClustIds = allWindOSPs.(dataMethod).(cellFieldName).clustIds;
            if ~isequal(savedClustIds, curClustIds)
                redo_current_cell = 1;
            end
        end
    end     
    
    
    nInit = switchh(gratingType, {'flashed', 'drifting'}, [3, 1]);
    blnk = cell(1, nInit);
    s_orig = struct;
    if haveSomeDataForThisCell && ~redo_current_cell
        s = allWindOSPs.(dataMethod).(cellFieldName);
        s_orig = s;
        for i = 1:length(dataCanCalculate)
            if ~isfield(s, dataCanCalculate{i})
                s.(dataCanCalculate{i}) = blnk;
            elseif length(s.(dataCanCalculate{i})) < nInit
                s.(dataCanCalculate{i}){nInit} = [];
            end                
        end
    else
        if exist('curClustIds', 'var')
            clustIdsFields = {'clustIds', {curClustIds}};
        else
            clustIdsFields = {};
        end
        structInit = [dataCanCalculate(:)'; repmat({{blnk}}, [1, length(dataCanCalculate)])  ];
        s = struct('l_bins', zeros(1, nInit, 'uint8'), 'r_bins', zeros(1, nInit, 'uint8'), 'windowProfile', {cell(1, nInit)}, ...
            structInit{:}, clustIdsFields{:} );
        assert(length(s) == 1);
%             'osps', {cell(1, nInit)},  'osps_odd', {cell(1, nInit)},  'osps_even', {cell(1, nInit)});        
    end
        
    if ~isfield(s, 'windowProfile')       
        s.windowProfile = cell(1, length(s.l_bins));
    end
    
    if ~isfield(s, 'lastUpdated')
        s.lastUpdated = zeros(1, length(s.l_bins));
    end
        
    idx_toRemove = find(strCcmp(dataToDelete, fieldnames(s)));
    if ~isempty(idx_toRemove)
%         error('!');
        %%
        for fld_i = 1:length(idx_toRemove)
            s = rmfield(s, dataToDelete{idx_toRemove(fld_i)});
        end                
    end    
        
    % reorder fields
    clustIdFld = iff(exist('curClustIds', 'var'), {'clustIds'}, {});
    firstFields = {'l_bins', 'r_bins', 'windowProfile'}; lastFields = {'STA', 'STA_oe', 'STA_jack', 'STA_odd_jack', 'STA_even_jack', 'meanRate', clustIdFld{:}, 'lastUpdated'}; %#ok<CCAT>
%     lastFields = setdiff(lastFields, dataToDelete); %%%%%% *******
    
    otherFields = setdiff(fieldnames(s), [firstFields, lastFields, dataToDelete]);
    fieldOrder = [firstFields, sort(otherFields(:)'), lastFields]; %#ok<TRSRT>
    
    s = orderfields(s, fieldOrder);        
    
%     [osp, osp_odd, osp_even]
    dataToCalc = dataToCalculate;
    
    w_idx = find( (s.l_bins == l_bin) & (s.r_bins == r_bin) & cellfun(@(c) isequal(c(:), windowProfile(:)), s.windowProfile), 1);
    
    dataIsOld = ~isempty(w_idx) && (redo_old_cell_data && s.lastUpdated(w_idx) < redo_cell_data_if_before);
    if ~isempty(w_idx) && ~redo_current_cell_data && ~dataIsOld
        isAvailable = cellfun(@(dn) isfield(s, dn) && (length(s.(dn)) >= w_idx) && ~isempty(s.(dn){w_idx}), dataToCalculate) ;

        dataToCalc = dataToCalc(~isAvailable);            
    
    elseif isempty(w_idx)
        w_idx = find(s.l_bins == 0, 1);
        if isempty(w_idx)
            w_idx = length(s.l_bins)+1;
        end
        
        s.l_bins(w_idx) = l_bin;
        s.r_bins(w_idx) = r_bin;       
        s.windowProfile{w_idx} = windowProfile;
        if exist('curClustIds', 'var')
            s.clustIds = curClustIds;
        end
    end
    
    
    if ~isempty(dataToCalc)
        osp_GC_useCorticalStateFromAllTrials = true;
%         fprintf('!');

        ospsToCalc = {};
        if any(strncmp('phaseTC_', dataToCalc, 8))
            if strcmp(responseType, 'raw')
                ospsToCalc = [ospsToCalc, 'osp_full'];
            elseif strcmp(responseType, 'gainCorrected')
                ospsToCalc = [ospsToCalc, 'osp_ph_oe', 'osp_ph_hoe'];
            end
        end
        
        if any(strCcmp({'stimF1oDCs', 'STA'}, dataToCalc)) % for these, we need 'osp_ph'
            ospsToCalc = [ospsToCalc, 'osp_ph']; 
        end

        if any(strcmp('STA_oe', dataToCalc)) % for these, we need 'osp_ph_oe'
            ospsToCalc = [ospsToCalc, 'osp_ph_oe']; 
        end
        
        if any(strCcmp({'tuningStats', 'F1oDC',   'STA_jack', 'STA_oe_jack', ...
                        'osp_ph_jackTrials', 'osp_ph_oe_jackTrials', 'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm'}, dataToCalc)) % for these, we need 'osp_full'
            ospsToCalc = [ospsToCalc, 'osp_full']; 
            if strcmp(responseType, 'gainCorrected') 
                ospsToCalc = [ospsToCalc, 'osp_full_oe'];
            end
        end
        
        ospsToCalc = unique([ospsToCalc(:)', dataToCalc(:)']);               
        
        osp_var_names = {'osp',  'osp_oe',    'osp_ph',  'osp_ph_oe',  'osp_ph_hoe',  'osp_full',   'osp_full_oe'};
        tfOspsToCalc = (strCcmp(osp_var_names, ospsToCalc)); 
        for i = find(tfOspsToCalc(:)')  % see which ones we already have.
            osp_name = osp_var_names{i};
            dontHaveThisOsp = ( (length(s.(osp_name)) < w_idx) || isempty(s.(osp_name){w_idx}) );
            tfOspsToCalc(i) = tfOspsToCalc(i) && (...
                   dontHaveThisOsp || ...
                  (redo_current_cell_data && any(strcmp(dataToCalc, osp_name)) ) || ...
                   redo_current_cell );
        end
        idxOspsToCalc = find(tfOspsToCalc); 
        nOspsToCalc = length(idxOspsToCalc);
        ospsToCalc = osp_var_names(idxOspsToCalc);
        
        
        haveMeanRate = (length(s.meanRate) >= w_idx) && ~isempty(s.meanRate{w_idx});
        if haveMeanRate
            meanRate = s.meanRate{w_idx};
        end
        if (isempty(allBins) || isempty(psthVals)) || ((~exist('meanRate', 'var') && any(strcmp('meanRate', dataToCalc))) && ~haveMeanRate);
            [allBins, psthVals, meanRate] = dbGetCellSpkStimHists(Gid, cellId);
            if ~haveMeanRate
                s.meanRate{w_idx} = meanRate;
            end
        end

        
        % if gain corrected - just calculated r_full_(raw) - from which
        % will calculate the other ones. 
        if nOspsToCalc > 0
            
            if strcmp(responseType, 'gainCorrected') 
                
                osp_full_nonGC = calcOspForPsthWindow(Gid, {allBins, psthVals}, l_bin, r_bin, false, windowProfile, 'osp_full', meanRate);
%                 if isempty(s.osp_full_nonGC{w_idx})                    
% %                     s.osp_full_nonGC{w_idx} = compress(single( osp_full_nonGC ));
%                 else
%                     osp_full_nonGC = decompress( s.osp_full_nonGC{w_idx} );
%                 end
                need_osp_all = osp_GC_useCorticalStateFromAllTrials && any(strCcmp({'osp_ph_oe', 'osp_ph_hoe', 'osp_full_oe'}, ospsToCalc));
                
                if (any(strCcmp({'osp_ph', 'osp_full'}, ospsToCalc)) || need_osp_all && isempty(s.osp_ph{w_idx})) || redo_current_cell
                    [osp_gainCorrected, osp_gainCorrected_full, stats] = getResponseCorrectedForCorticalState(Gid, cellId, osp_full_nonGC, 'all');
                    stats.R = compress(single(osp_gainCorrected));
                    s.osp_ph{w_idx} = stats;
                    s.osp_full{w_idx} = compress(single(osp_gainCorrected_full));
                end
                
                if need_osp_all
                    stats = s.osp_ph{w_idx};
                    corticalGainParams = stats.corticalGainParams;
                else
                    corticalGainParams = [];
                end
                
                if any(strcmp({'osp_ph_oe', 'osp_full_oe'}, ospsToCalc)) || redo_current_cell
                    [osp_gainCorrected_oe, osp_gainCorrected_full_oe, stats_oe] = getResponseCorrectedForCorticalState(Gid, cellId, osp_full_nonGC, 'oe', corticalGainParams); %#ok<NASGU>
%                     stats_oe.R = compress(single( osp_gainCorrected_oe ));
%                     s.osp_ph_oe{w_idx} = stats_oe;
                    s.osp_ph_oe{w_idx} = compress(single( osp_gainCorrected_oe ));
                    s.osp_full_oe{w_idx} = compress(single( osp_gainCorrected_full_oe ));
                    
                end
                if any(strcmp({'osp_ph_hoe'}, ospsToCalc)) || redo_current_cell
                    [osp_gainCorrected_hoe, ~, stats_hoe] = getResponseCorrectedForCorticalState(Gid, cellId, osp_full_nonGC, 'hoe', corticalGainParams); %#ok<NASGU>
%                     stats_hoe.R = compress(single( osp_gainCorrected_hoe ));
%                     s.osp_ph_hoe{w_idx} = stats_hoe;
                    s.osp_ph_hoe{w_idx} = compress(single( osp_gainCorrected_hoe ));
                end
                    

            elseif strcmp(responseType, 'raw') 
        
                [osp_vars{1:nOspsToCalc}] = calcOspForPsthWindow(Gid, {allBins, psthVals}, l_bin, r_bin, false, windowProfile, osp_var_names(idxOspsToCalc), meanRate);
    %             osp_full_tmp = calcOspForPsthWindow(Gid, {allBins, psthVals}, l_bin, r_bin, false, windowProfile, {'osp_full'}, meanRate);

                for i = 1:nOspsToCalc
                    s.(osp_var_names{idxOspsToCalc(i)}){w_idx} = compress( single(osp_vars{i}) );
                end
            end
        end
        
        if strcmp(responseType, 'raw')
            get_osp_fld         = @(fld) double( decompress( s.(fld){w_idx}) );
            get_osp_fld_oe_full = @(fld) double( decompress( s.(fld){w_idx}) );
        else
            get_osp_fld         = @(fld) double( decompress( s.(fld){w_idx}.R) );
            get_osp_fld_oe_full = @(fld) double( decompress( s.(fld){w_idx}) );
        end

        
        if any(strncmp('phaseTC_', dataToCalc, 8)) % ... || any(strCcmp({'phaseTC_CCs', 'phaseTC_CC_ps', 'phaseTC_Dots'}, dataToCalc))
            if  strcmp(responseType, 'raw') 
                R_full = decompress (s.osp_full{w_idx});
                nTrials = size(R_full,4);
                idx_odd = odd([1:nTrials]);  idx_even = ~idx_odd;
                R_odd = mean(R_full(:,:,:,idx_odd), 4);
                R_even = mean(R_full(:,:,:,idx_even), 4);
                
                idx_halfOdd = odd(floor([1:nTrials] / 2)); idx_halfEven = ~idx_halfOdd;
                R_halfOdd = mean(R_full(:,:,:,idx_halfOdd), 4);
                R_halfEven = mean(R_full(:,:,:,idx_halfEven), 4);
                
                idx_firstHalf = [1:nTrials] <= nTrials/2; idx_secondHalf = ~idx_firstHalf;
                R_firstHalf = mean(R_full(:,:,:,idx_firstHalf), 4);
                R_secondHalf = mean(R_full(:,:,:,idx_secondHalf), 4);
            else
%                 R_oe = decompress (s.osp_ph_oe{w_idx}.R);
                R_oe = decompress (s.osp_ph_oe{w_idx});
                R_odd = R_oe(:,:,:,1);
                R_even = R_oe(:,:,:,2);
                
                R_hoe = decompress (s.osp_ph_hoe{w_idx});
                R_halfOdd = R_hoe(:,:,:,1);
                R_halfEven = R_hoe(:,:,:,2);
                
                % skip first/second half for now
            end
                
            
            if any(strcmp('phaseTC_CCs_oe', dataToCalc)) || any(strcmp('phaseTC_CC_ps_oe', dataToCalc))
                [phaseTC_CCs_oe, phaseTC_CC_ps_oe] = getAllPhaseTcCCs(R_odd, R_even);                
                s.phaseTC_CCs_oe{w_idx} = single(phaseTC_CCs_oe);
                s.phaseTC_CC_ps_oe{w_idx} = single(phaseTC_CC_ps_oe);
            end
            if any(strcmp('phaseTC_CCs_hoe', dataToCalc)) || any(strcmp('phaseTC_CC_ps_hoe', dataToCalc))
                [phaseTC_CCs_hoe, phaseTC_CC_ps_hoe] = getAllPhaseTcCCs(R_halfOdd, R_halfEven);
                s.phaseTC_CCs_hoe{w_idx} = single(phaseTC_CCs_hoe);
                s.phaseTC_CC_ps_hoe{w_idx} = single(phaseTC_CC_ps_hoe);
            end
            if any(strcmp('phaseTC_CCs_fs', dataToCalc)) || any(strcmp('phaseTC_CC_ps_fs', dataToCalc))
                [phaseTC_CCs_fs, phaseTC_CC_ps_fs] = getAllPhaseTcCCs(R_firstHalf, R_secondHalf);
                s.phaseTC_CCs_fs{w_idx} = single(phaseTC_CCs_fs);
                s.phaseTC_CC_ps_fs{w_idx} = single(phaseTC_CC_ps_fs);
            end

%             if any(strcmp('phaseTC_Dots', dataToCalc))
%                 phaseTC_Dots = getAllPhaseTcCCs_sampled(stimPSTH_vals_allTrials, l_bin-binOffset, r_bin-binOffset, 'dot');                            
%                 s.phaseTC_Dots{w_idx} = single(phaseTC_Dots);
%             end                        
        end
                
        if any(strCcmp({'stimF1oDCs', 'stimF1oDC_jackStds'}, dataToCalc))
%             r_ph = get_osp_fld('osp_ph');
%             r_ph = get_osp_fld('osp_ph');
%             [nOri, nSp, ~] = getGratingStimType(Gid);       
%             r_ph_rows = permuteReshapeData(r_ph, 3)';
%             stimF1oDCs = getF1oDC(r_ph_rows);
%             stimF1oDCs = reshape(stimF1oDCs, [nOri, nSp]);                    
            r_full = get_osp_fld_oe_full('osp_full');
            [stimF1oDCs, stimF1oDC_jackStds] = getStimF1oDCs(r_full);

            s.stimF1oDCs{w_idx} = single(stimF1oDCs);
            s.stimF1oDC_jackStds{w_idx} = single(stimF1oDC_jackStds);
            
        end

        if any(strcmp('STA', dataToCalc))                                    
            r_ph = get_osp_fld('osp_ph');                
            STA = getSTAfromOSP(Gid, r_ph);            
            s.STA{w_idx} = single(STA);                        
        end
        
        if any(strcmp('STA_oe', dataToCalc))
            r_ph_oe = get_osp_fld_oe_full('osp_ph_oe');
            STA_odd = getSTAfromOSP(Gid, r_ph_oe(:,:,:,1));
            STA_even = getSTAfromOSP(Gid, r_ph_oe(:,:,:,2));
            s.STA_oe{w_idx} = cat(3, single(STA_odd), single(STA_even));                        
        end
        
        if any(strCcmp({'STA_jack', 'STA_odd_jack', 'STA_even_jack'}, dataToCalc))
            stimInfo = getGratingStimType(Gid);
            nJackSegments_STA = switchh(stimInfo.nTrials, [4, 16,   10], [4, 4,  5]);
            jackknifeMethod_STA = 'trials';
        
            if strcmp(responseType, 'raw')
                r_full_all = get_osp_fld_oe_full('osp_full');
                r_full_oe  = get_osp_fld_oe_full('osp_full');
            elseif strcmp(responseType, 'gainCorrected')
                r_full_all = get_osp_fld_oe_full('osp_full');
                r_full_oe  = get_osp_fld_oe_full('osp_full_oe');               
            end
            
%             r_full = double( decompress(  s.osp_full{w_idx}) );
            if any(strcmp('STA_jack', dataToCalc))
                s.STA_jack{w_idx} = getJackknifedSTA(r_full_all, Gid, 'all', nJackSegments_STA, jackknifeMethod_STA);
            end

            if any(strcmp('STA_odd_jack', dataToCalc))
                s.STA_odd_jack{w_idx} = getJackknifedSTA(r_full_oe, Gid, 'odd', nJackSegments_STA, jackknifeMethod_STA);
            end

            if any(strcmp('STA_even_jack', dataToCalc))
                s.STA_even_jack{w_idx} = getJackknifedSTA(r_full_oe, Gid, 'even', nJackSegments_STA, jackknifeMethod_STA);
            end
        end
                
        if any(strCcmp({'osp_ph_jackTrials', 'osp_ph_oe_jackTrials', 'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm'}, dataToCalc))
%             r_full = single ( decompress(  s.osp_full{w_idx}) );
    
            if strcmp(responseType, 'raw')
                r_full_all = get_osp_fld_oe_full('osp_full');
                r_full_oe  = get_osp_fld_oe_full('osp_full');
            elseif strcmp(responseType, 'gainCorrected')
                r_full_all = get_osp_fld_oe_full('osp_full');
                r_full_oe  = get_osp_fld_oe_full('osp_full_oe');               
            end
            
            doSmoothingNow = any(strCcmp({'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm'}, dataToCalc));

            if doSmoothingNow
                [smoothPhase_Width, smoothPhase_Method] = curPhaseSmoothing; 
                assert( (flashedOrDrifting(Gid) == 2) && strcmp(curCmpType(''), 'phase') && strcmp(curGratingType(''), 'drifting') );
                smoothPhasesFunc = @(R) smoothOSP_Phs(R, smoothPhase_Method, smoothPhase_Width);
            end
            %%
            if any(strCcmp({'osp_ph_jackTrials', 'osp_ph_jackTrials_sm'}, dataToCalc))
                R_jackTrials_aa = getPhaseTuningJackknifedTrials(r_full_all, 'aa');
                if any(strcmp({'osp_ph_jackTrials_sm'}, dataToCalc)) && doSmoothingNow
                    R_jackTrials_aa_smoothed = cellfun(smoothPhasesFunc, R_jackTrials_aa, 'un', 0); % smooth
                    s.osp_ph_jackTrials_sm{w_idx} = compress(R_jackTrials_aa_smoothed);
                end
                s.osp_ph_jackTrials{w_idx} = compress(R_jackTrials_aa);
            end
%                 R_jackTrials_oe_d_sm = cellfun(smoothPhasesFunc, R_jackTrials_oe_d, 'un', 0); % smooth    
            if any(strCcmp({'osp_ph_oe_jackTrials', 'osp_ph_oe_jackTrials_sm'}, dataToCalc))
                R_jackTrials_oe = getPhaseTuningJackknifedTrials(r_full_oe, 'oe');
                if any(strcmp({'osp_ph_oe_jackTrials_sm'}, dataToCalc)) && doSmoothingNow
                    R_jackTrials_oe_smoothed = cellfun(smoothPhasesFunc, R_jackTrials_oe, 'un', 0); % smooth
                    s.osp_ph_oe_jackTrials_sm{w_idx} = compress(R_jackTrials_oe_smoothed);
                end
                s.osp_ph_oe_jackTrials{w_idx} = compress(R_jackTrials_oe);
            end
            
            
            
        end
                
        if any(strcmp('osp_pred', dataToCalc))                                                
%             STA = double(s.STA{w_idx});
            MID_filename = mid_getPreferredMIDfile(Gid, cellId); % getName('MID_file', 4470, 2, downSampFactor, frameMode);            
            if exist(MID_filename, 'file')
                S = load(MID_filename);                
            else
                error('MID unavailable');
            end            
            r_ph = get_osp_fld('osp_ph');
            
            [R_pred, nlin_func, cc] = predictOSPfromMID(Gid, S.MID, r_ph);
            s.osp_pred{w_idx} = {R_pred, nlin_func, cc};
        end
                

        if any(strcmp('tuningStats', dataToCalc))            
            sd = siteDataFor(Gid);
            if strncmp(sd.stimType, 'Movie:Flashed_Gratings', 12)
                nSamplesAv = r_bin-l_bin+1;
            else % drifting gratings
                bckgBinSize_ms = 4+1/6;
                utf_Hz = 1/sd.tempPeriod_sec;
                nSamplesAv = round((utf_Hz*1000)/bckgBinSize_ms);        
            end
            bckgSkip_ms = 250;
            bckgSamples = getBackgroundSpikes(Gid, cellId, bckgSkip_ms, nSamplesAv);
            
            osp_full = get_osp_fld_oe_full('osp_full');

            if ~exist('meanRate', 'var') 
                meanRate = getOspDataForPsthWindow(Gid, cellId, [], [], l_bin, r_bin, windowProfile, {'meanRate'});
            end
            osp_full = osp_full / mean(osp_full(:)) * meanRate;
            
            s.tuningStats{w_idx} = calcDegreeOfTuningStats(osp_full, bckgSamples, Gid, cellId);
        end
        
        if any(strcmp('windowStats', dataToCalc))
            allWindowStats = getPSTHwindowData(Gid, cellId);
            allWindowStats_curWindow = structfun(@(X) X(r_bin, l_bin), rmfield(allWindowStats, 'clustIds'), 'un', 0);
            s.windowStats{w_idx} = allWindowStats_curWindow;
        end
                
        if ~curMatchDB && strcmp(cellClustStr, 'cells')
            curClustIds = getClusterCellAssignment(Gid, cellId);
            s.clustIds = curClustIds;
        end
    end
    
    if ~isequaln(s, s_orig)
        %%
        s.lastUpdated(w_idx) = now;
        allWindOSPs.(dataMethod).(cellFieldName) = s;
        saveCount = saveCount+1;
        saveFile = true;
    end

    
%     curClustIds = getClusterCellAssignment(Gid, cellId);
%     allWindOSPs.(dataMethod).(cellFieldName).clustIds = curClustIds;
    
    dataKeepCompressed = {'osp_ph_jackTrials', 'osp_ph_oe_jackTrials', 'osp_ph_jackTrials_sm', 'osp_ph_oe_jackTrials_sm'};  % don't decompress these

    % OUTPUT the requested data
    for i = 1:length(dataToCalculate)
        argout_i = s.(dataToCalculate{i}){w_idx};
        if isstruct(argout_i)  && ~any(strcmp(dataToCalculate{i}, dataKeepCompressed)) 
            argout_i = decompress(argout_i);
        end
        varargout{i} = argout_i;
    end
    
    

    if saveFile && (saveCount > saveCountSpacing)
        %%
        allCellData  = orderfields(allWindOSPs.(dataMethod));         %#ok<NASGU>
        fprintf('\n[saving osp data (%s)... ', dataMethod); tic;
        save(ospWindowDataFile, '-struct', 'allCellData', '-v6');    
        t_sec = toc; fprintf('done (%.2f s)\n', t_sec);        
        saveCount = 0;                
    end                        
                
    
end
    


% function s = cellFieldName_fcn(Gid, cellId, methName) 
%     s = sprintf('Gid_%04d_cell_%d_%s', Gid, cellId, methName);            
% end


%{
            [bins_allTrials, stimPSTH_vals_allTrials] = dbGetCellSpkStimHists(Gid, cellId, struct('trialGrouping', 'individual'));
            binOffset = find( abs(allBins - bins_allTrials(1)) < 1e-5, 1) - 1;
            nTrials = size(stimPSTH_vals_allTrials,3);            
            stimPSTH_vals_allTrials = reshape(stimPSTH_vals_allTrials, [length(bins_allTrials), nOri, nSp, nPh, nTrials]);
%             phaseTC_cc_vars = strCcmp({'phaseTC_CCs', 'phaseTC_CC_ps'}, dataToCalc);
            
            if any(strcmp('phaseTC_CCs_oe', dataToCalc)) || any(strcmp('phaseTC_CC_ps_oe', dataToCalc))
                [phaseTC_CCs_oe, phaseTC_CC_ps_oe] = getAllPhaseTcCCs_sampled(stimPSTH_vals_allTrials, l_bin-binOffset, r_bin-binOffset, 'cc', 'odd_vs_even');                
                s.phaseTC_CCs_oe{w_idx} = single(phaseTC_CCs_oe);
                s.phaseTC_CC_ps_oe{w_idx} = single(phaseTC_CC_ps_oe);
            end

%}