function varargout = getPSTHwindowData(Gid, cellId, statTypes, psthWindow_ms, psthMethod, l_bins, r_bins)

    %%%%% parameters %%%%
    global psthStatsSettings
    persistent allPsthWindowData saveCount;

    if isempty(psthStatsSettings)
        setPsthGlobalSettings;
    end

    checkWhetherCellAssignmentHasChanged = false;
    
    redoAllCells = false;
    redoCurrentCell = false;
    redoCurrentCellStats = false;
    ignoreFile = false;
    
    saveCountSpacing = 30;
    
    decompressOnInitialLoad = true;
    assumeWeHaveAllStats = exist('Gid', 'var') && length(Gid) > 1;
    removeUnnecessaryStats = false;    
    %%%%%%%%%%%%%%%%%%%%%
    
    stg = psthStatsSettings;
    assert( ~(stg.compressToList && stg.compressToSparse) ); % maximum of 1 of these should be true;
    assert( ~stg.compressToSparse || (stg.compressToSparse && strcmp(stg.storeUnassignedValuesAs, 'zeros') && stg.storeLogOfPvalues ) );     % if compressing to sparse then should store unassigned as zeros, and therefore should take log of pvals.    
    assert( ~stg.storeLogOfPvalues || (stg.storeLogOfPvalues && (stg.nudge0pval_amt > 0 || strcmp(stg.storeUnassignedValuesAs, 'nan')) )  );  % if storeLogOfPvalues = false, should storeUnassigned as nans.        
    assert( ~stg.compressToSparse || (stg.compressToSparse && strcmp(stg.statsPrecision, 'double')) );  % if storing as sparse, must be double
     
    if ~exist('Gid', 'var') || isempty(Gid)
        gratingType = curGratingType('');
    elseif ischar(Gid) && strcmp(Gid, 'save')
        if isempty(allPsthWindowData) % || (saveCount == 0)
            return;
        end
        fn = fieldnames(allPsthWindowData);
        isFlashed  = ~isempty(strfind(fn{1}, '_f_'));
        isDrifting = ~isempty(strfind(fn{1}, '_d_'));
        assert(xor(isFlashed, isDrifting));
        gratingType = iff(isFlashed, 'flashed', 'drifting');        
    elseif isnumeric(Gid)                
        gratingType = flashedOrDrifting(Gid, 'str');  
    end    
    
    if ~exist('psthMethod', 'var') || isempty(psthMethod)
        psthMethod = 'stimulus';
    end
    psthMethod = psthMethod(1:2);
    ospFcn_str = stg.ospPhCompressFcn;        
    
%     redoCurrentCellStats_if_short = false;
    if isempty(saveCount)
        saveCount = 0;
    end

    noGidInput = ~exist('Gid', 'var') || isempty(Gid);        
    
    if ~exist('statTypes', 'var') || isempty(statTypes)
%         statTypes = {'rep_p', 'rep_tstat', 'rho', 'rho_t', 'rho_p',  'r_var', 'r_entropy'};
%         statTypes = {'r_var'};
%         statTypes = {'rep_p', 'rep_tstat', 'rho', 'rho_p', 'rho_t',  'tau', 'tau_p', 'tau_t', 'r_var', 'r_entropy'};
%         statTypes = {'rep_p', 'rep_tstat',  'rho', 'rho_p', 'rho_p_nz', 'rho_p_nznz', 'r_entropy', 'n_unique'};
%         statTypes = {'rep_p', 'rep_tstat',  'cc_p',  'tau_p', 'tau_t',  'rho_p', 'rho_p_nz', 'rho_p_nznz', 'r_entropy'};
        statTypes = {'cc_p'};
    end
    if ischar(statTypes)
        statTypes = {statTypes};
    end
    
    if ~exist('psthWindow_ms', 'var') || isempty(psthWindow_ms)
        switch gratingType
            case 'flashed',  psthWindow_ms = [-300, 200]; 
            case 'drifting', psthWindow_ms = [0, 8+1/3];
        end    
    end
        
        
    if nargin < 6
        l_bins = [];
        r_bins = [];
    end
    
    cellClustStr = curGroupingType('');    
    fet_str = iff(curMatchDB, '_DB', ['_' curSortingFeatures('')]);
            
    psthType = [cellClustStr fet_str '_' gratingType(1) '_' psthMethod '_' ospFcn_str ];


    
    psthWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'allPsthWindowData_' psthType '.mat'];
    
    if ischar(Gid) && strcmp(Gid, 'save') 
        if ~isempty(allPsthWindowData.(psthType))
            %%
            fprintf('[saving Psth-Window Data (%d new entries) (%s) ...', saveCount, psthType); tic;
            tmp = allPsthWindowData.(psthType); %#ok<NASGU>
            save(psthWindowDataFile, '-struct', 'tmp', '-v6');    
            t_sec = toc; fprintf('done (%.2f s)]\n', t_sec);        
            saveCount = 0;            
        end
        return;
    end

    if isempty(allPsthWindowData)
        allPsthWindowData = struct;
    end
        
    if ~isfield(allPsthWindowData, psthType) || redoAllCells
        
        if ~exist(psthWindowDataFile, 'file') || ignoreFile || redoAllCells 
            allPsthWindowData.(psthType) = struct;
        else                
            method_str = iff(strcmp(psthMethod, 'st'), 'stimulus', 'spike');
            fprintf('Loading %s %s grating Psth-Window Data using %s method, taking %s over phases.  ... ', cellClustStr, upper(gratingType), upper(method_str), upper(ospFcn_str) ); tic;             
            
            allPsthWindowData.(psthType) = load(psthWindowDataFile);
%             if decompressOnInitialLoad
%                 decompressStats = @(s) structfun(@decompressStat, s, 'un', 0);
%                 allPsthWindowData = structfun(decompressStats, allPsthWindowData, 'un', 0);
%             end               
            fprintf(' done. '); toc;
        end
    end

    
    if noGidInput  % do for all gids/cellids in flashed grating cells file
        sortByRep = false;
%         S1 = load('cellsGroups_movie_fg.mat');
%         frm = [S1.movieGroups_fg.frameLength_ms];
%         grp_idx = find ( frm > 90 );
%         all100Gids = [S1.movieGroups_fg(grp_idx).Gid];
% 
%         S = load('flashedGratingCells_all.mat');        
%         idx2 = find( arrayfun(@(g_id) any(all100Gids == g_id), [S.allCells.Gid]) );        
%         Gid = [S.allCells(idx2).Gid];
%         cellId = [S.allCells(idx2).cellId];            

%         S1 = load('cellsGroups_movie_fg.mat');
%         st = {S1.movieGroups_fg.stimType};
%         grp_idx = cellfun(@(s) strcmp(s, 'Movie:Flashed_Gratings:36x10x8(16x1)'), st );
%         all2x8Gids = [S1.movieGroups_fg(grp_idx).Gid];            
% 
%         S = load('flashedGratingCells_all.mat');        
%         idx2 = find( arrayfun(@(g_id) any(all2x8Gids == g_id), [S.allCells.Gid]) );        
%         Gid = [S.allCells(idx2).Gid];
%         cellId = [S.allCells(idx2).cellId];            

        S = load('flashedGratingCells_all.mat'); allCells = S.allCells;
        if sortByRep
            rep = arrayfun(@(s) s.stats.rep_ori_sp_maxPhase_str, allCells);
            allCells = allCells(ord(rep, 'descend'));
        end    
        
        Gid = [allCells.Gid];
        cellId = [allCells.cellId];              
%         Gid = [4470 4470];
%         cellId = [0 2];
%         Gid = [4494, 4476, 4522, 4546, 5112];
%         cellId = [0, 0, 3, 0, 3];

%         Gid = [4462, 4466, 4476, 4476, 4508, 4706, 4718, 5112];
%         cellId = [1 3, 4, 5, 0, 1, 0, 0];
        
%         Gid = [4470, 4470, 4476, 4482, 4488, 4506, 4508, 4706, 4712, 4796];
%         cellId = [1, 4, 1, 4, 0, 1, 1, 0, 1, 0        ];
        
        
%         S = load('flashedGratingCells_all.mat');
%         allCells = S.allCells;
%         
%         ignoreMUs = false;
%         if ignoreMUs
%             allCells = allCells( [allCells.cellId] > 0);
%         end
%         Gid = [allCells.Gid];
%         cellId = [allCells.cellId];                            
    end    

    % if doesn't exist in loaded file, calculate it
    if assumeWeHaveAllStats && length(Gid) > 10
        cellnames = arrayfun(@(g,c) cellFieldName_fcn(g, c, psthMethod), Gid(:), cellId(:), 'un', 0);
        fn = fieldnames(allPsthWindowData.(psthType));
        data_C = struct2cell(allPsthWindowData.(psthType));
        if ~isequal(fn, cellnames)
            
            if length(fn) == length(cellnames)
                allPsthWindowData.(psthType) = orderfields(allPsthWindowData.(psthType), cellnames);
            elseif length(fn) > length(cellnames)

                idxs = idxsOfBinA(fn, cellnames);
                data_C = data_C(idxs);
            end
            3;
        end        
        
%         fn = fieldnames(data_C{1});
%         fn_remove = setdiff(fn, statTypes);        
%         data_C = cellfun(@(c) rmfield(c, fn_remove), data_C, 'un', 0);
        
        varargout = {data_C};
        return;
    end
    
    
    if length(Gid) > 1
        progressBar('init-', length(Gid), 60);
    end
    saveFile = false;
    for cell_i = 1:length(Gid)  
        if length(Gid) > 1
            progressBar;
        end
        psth_name = cellFieldName_fcn(Gid(cell_i), cellId(cell_i), psthMethod );
        
        haveSomeDataForThisCell = isfield(allPsthWindowData.(psthType), psth_name) ;
        
        if haveSomeDataForThisCell && ~curMatchDB && strcmp(cellClustStr, 'cells') && checkWhetherCellAssignmentHasChanged 
            curClustIds = getClusterCellAssignment(Gid, cellId);
            savedClustIds = allPsthWindowData.(psthType).(psth_name).clustIds;
            if ~isequal(savedClustIds, curClustIds)
                redoCurrentCell = 1;
            end
        end        
        
        if haveSomeDataForThisCell && ~redoCurrentCell
            s = allPsthWindowData.(psthType).(psth_name); 
        else
            s = [];
        end
        
        %{ 
        if redoCurrentCellStats_if_short
            if haveSomeDataForThisCell        %&& redoCurrentCellStats_if_short
                fn = fieldnames(s);
                tooShort = (length( s.(fn{1}) ) < 3690 ); % all windows <= 200 ms for [-300, 200] window
            else
                tooShort = false;
            end
        end
        %}
    
%         s_startWith = iff(isempty(l_bins), [], s);
%         s_startWith = s;%
        
        if assumeWeHaveAllStats
            statToDo = {};
        elseif haveSomeDataForThisCell && ~redoCurrentCell && ~redoCurrentCellStats && isempty(l_bins) % && ~tooShort
            if ~isempty(s)
                statToDo = setdiff(statTypes, fieldnames(s));
            end
        else
            statToDo = statTypes;            
        end
        
        if ~isempty(statToDo)            
            D = getStatsData(Gid(cell_i), cellId(cell_i), statToDo, psthWindow_ms, ~noGidInput, l_bins, r_bins, s);
            for stat_i = 1:length(statToDo)
                s.(statToDo{stat_i}) = D{stat_i};
            end
            if noGidInput
                                
            end
            
            if ~curMatchDB && strcmp(cellClustStr, 'cells')
                curClustIds = getClusterCellAssignment(Gid, cellId);
                s.clustIds = curClustIds;
            end

            allPsthWindowData.(psthType).(psth_name) = s;            
            saveFile = true;
            saveCount = saveCount +1;
%             fprintf('*')
        end
%             [tmp1, tmp2, areas] = getStatsData(Gid(cell_i), cellId(cell_i));
%             allPsthWindowData_S.(psth_name) = struct('pvals', s.pvals, 'tstats', s.tstats, 'areas', areas);
%             saveFile = true;            
    end
    if length(Gid) > 1
        progressBar('done')
    end
    if saveFile && ~ignoreFile && (saveCount > saveCountSpacing)        
        tmp  = orderfields(allPsthWindowData.(psthType));         %#ok<NASGU>
        fprintf('saving psth window data ... '); tic;
        save(psthWindowDataFile, '-struct', 'tmp', '-v6');    
        t_sec = toc; fprintf('done (%.2f s)\n', t_sec);        
        saveCount = 0;
    end                        
    
            
    if length(Gid) == 1
        if ~exist('s', 'var')
            psth_name = cellFieldName_fcn(Gid, cellId, psthMethod);
            s = allPsthWindowData.(psthType).(psth_name);
        end
        
%         curClustIds = getClusterCellAssignment(Gid, cellId);
%         allPsthWindowData.(psthType).(psth_name).clustIds = curClustIds;

        
        if removeUnnecessaryStats
            flds_dontWant = setdiff(fieldnames(s), statTypes);
            if ~isempty(flds_dontWant); % ie only want subset of the fields
                s = rmfield(s, flds_dontWant);            
            end
        end
        flds = fieldnames(s);    
        if ((stg.compressToList && isvector( s.(flds{1}) )) || stg.compressToSparse) && ~decompressOnInitialLoad 
            s = structfun(@decompressStat, s, 'un', 0);            
        end
        s = structfun(@numericFull, s, 'un', 0);

        if nargout <= 1
            varargout = {s};
        elseif nargout == length(flds)
            varargout = struct2cell(s); 
        else
            error('Incorrect # of outputs');
        end
    else
        varargout = {};
    end
        
    
end
    
function y = numericFull(x)
    if isnumeric(x) && issparse(x)
        y = full(x);
    else
        y = x;
    end

end



function stats = getStatsData(Gid, cellId, statTypes, psthWindow_ms, doShort, l_bins, r_bins, stats_S)

    global psthStatsSettings
    stg = psthStatsSettings;        

    nStats = length(statTypes);
    statTypesToDo = statTypes;
    idx_stats = 1:nStats;
    weightByPsth = false;
%     idx_area = find(strcmp('area', statTypes),1);
%     doArea = ~isempty(idx_area);
%     if doArea
%         statTypesToDo(idx_area) = [];
%     end
%     idx_stats = setdiff(1:nStats, idx_area);
    
    psthMethod = 'stimulus';
    trialGrouping = 'odd/even';
    splitCphTrials = true;
        
    histOpts = struct('splitWindowIfCph', splitCphTrials, 'psthMethod', psthMethod, 'trialGrouping', trialGrouping, 'psthWindow_ms', psthWindow_ms);
    [PSTH_bins, stimPSTH_vals_oe] = dbGetCellSpkStimHists(Gid, cellId, histOpts);    
    
    initFcn = str2func(stg.storeUnassignedValuesAs);
    if flashedOrDrifting(Gid) == 1
        binWidth = 25/6;
        saveOsp = false;
    else
        binWidth = 25/3; %***
        saveOsp = true;
    end
    minWindowWidth_bins = 1;
    if doShort
        maxWindowWidth_ms = 120;
%         maxWindowWidth_ms = 35;
    else
        maxWindowWidth_ms = 35;
%         maxWindowWidth_ms = 200;
    end
    maxWindowWidth_bins = ceil( maxWindowWidth_ms/binWidth );
    
    nBins = length(PSTH_bins);
%     binStart = 1; %find(PSTH_bins > 15,1, 'first');
%     binEnd   = nBins ;%find(PSTH_bins < 120,1, 'last');

    binStart = find(PSTH_bins > 0,1, 'first');
    binEnd   = nBins ;%find(PSTH_bins < 120,1, 'last');
    
    
%     [uori, usp, uph] = dbGetUniqueOriSpPh('Gid', Gid);
%     [nOri, nSp, nPh] = deal(length(uori), length(usp), length(uph));

       
%     if NStim == 0

    % guess: take top 20 stimuli for PSTH (ordered by mean .* max of first 120 ms)        
    if weightByPsth
        stimPSTH_vals = mean(stimPSTH_vals_oe,3); 
        idx_relevant = (PSTH_bins > 20) & (PSTH_bins < 120);
    %     idx_relevant = (PSTH_bins < 120);                            
        psth_means = mean(stimPSTH_vals(idx_relevant,:), 1);
        psth_maxes = max(stimPSTH_vals(idx_relevant,:), [], 1);                
        new_order = ord(psth_means .* psth_maxes, 'descend');
        stimPSTH_vals_ordered = stimPSTH_vals(:,new_order);
        nStimMax = 20;        
        curPsth = mean(stimPSTH_vals_ordered(:,1:nStimMax), 2);        
        curPsth = rectified( curPsth - mean(curPsth));       
    else
        curPsth = [];
    end
%         stimPSTH_vals_cum = cummean(stimPSTH_vals_ordered,2);
%         stimPSTH_sums = sum( stimPSTH_vals_ordered, 1);
    %find(stimPSTH_sums / stimPSTH_sums(1) <= nStimTh, 1);

%         [pvals, tstats, areas]  = deal( nan(nBins, nBins, precn) );

    precn = psthStatsSettings.statsPrecision;
    
    stats = cell(1,nStats);
    if ~isempty(stats_S) 
        fn = fieldnames(stats_S);
        for i = 1:length(statTypes)
            idx_stat = find(strcmp(statTypes{i}, fn), 1);  % % 'un'=0 in case some are empty.
            if isempty(idx_stat)
                stats{i} = initFcn(nBins, nBins, precn);
            else
                stats{i} = stats_S.(statTypes{i}); 
            end
        end
    else        
        stats(:) = {initFcn(nBins, nBins, precn)};        % initialize as double for double precision.
    end
    
%     stimType = getGratingStimType(Gid);
    shiftCenter = true;        
        
    CellId = [];
    if isempty(l_bins)  %% do all bins of a
        n_approx = length(binStart:binEnd)*minWindowWidth_bins;    
        l_bins = zeros(1, n_approx); r_bins = zeros(1, n_approx);
        idx = 1;
        for l_bin = binStart : binEnd        
            r_binEnd = min(binEnd, l_bin+maxWindowWidth_bins -1);
            for r_bin = l_bin + (minWindowWidth_bins-1) : r_binEnd                                            
                l_bins(idx) = l_bin;
                r_bins(idx) = r_bin;
                idx = idx + 1;
            end
        end
        l_bins(idx:end) = []; r_bins(idx:end) = [];
        LR_unique = unique([l_bins(:), r_bins(:)], 'rows');
        if length(l_bins) ~= size(LR_unique, 1);
            error('not unique ');
        end
        
    else
        
        %GidCellId = [Gid, cellId];
    end
    
    if saveOsp
        CellId = cellId;
    end
    
%     idx = sub2indV( length(PSTH_bins)*[1, 1], [l_bins(:), r_bins(:)]);
%     idx_dontHave
    
%     assert(length(l_bins) == idx-1);
        
    for wind_i = 1:length(l_bins);    
        L = l_bins(wind_i); R = r_bins(wind_i);
        stat_tmp = getOspStatVsPsthBinning(Gid, PSTH_bins, stimPSTH_vals_oe, L, R, shiftCenter, curPsth, statTypesToDo, CellId);
        for stat_i = 1:length(stat_tmp)
            stats{idx_stats(stat_i)}(R, L) = stat_tmp(stat_i);
        end                
    end
    
    if stg.storeLogOfPvalues
        pval_stats_idx = find( cellfun(@(s) ~isempty(strfind(s, '_p')), statTypesToDo) );
        for i = pval_stats_idx(:)'
%             stats{i} = pval2NegLogPval( stats{i} );
           stats{i} = fixNegLogPval( stats{i}, psthStatsSettings.nudge0pval_amt );
        end
    end
    
    if stg.compressToList
        stats = cellfun(@compressOutNans, stats, 'un', 0);
    elseif stg.compressToSparse    
        if ~issparse(stats{1})
            stats = cellfun(@sparse, stats, 'un', 0);
        end
    end
    
    
end

function X = decompressStat(X)
    needToDecompressIfSparse = false;

    if isvector(X)
        X = decompressOutNans(X);
    elseif issparse(X) && needToDecompressIfSparse
        X = full(X);
    end
end
        
function s = cellFieldName_fcn(Gid, cellId, psthMethod) 
    s = sprintf('Gid_%04d_cell_%d_%s', Gid, cellId, psthMethod);
end


function p = fixNegLogPval(p, amt)    
    % we assigned p==0 a small amount of 1e-100. But this is a huge outlier
    % if the other p values are not so significant. so now we set those
    % pvalues to more intermediate values (close to the original)
    defaultAmt_idx = (p == -log10(amt) );
    if all(p) % were no values other than p == 0.
        return;
    end
    realMaxP = max(  p(~defaultAmt_idx(:))  );
    p(defaultAmt_idx) = ceil(realMaxP) + 1;
    
    % we also assigned 'logp'-values of -1 to p-values of 1. This can be too
    % negative- we want it to be only slightly below 0 compared with the
    % amounts above zero.
    r = .01;
    maxP = max(p(p ~= 0));
    if maxP < 0, maxP = -1; end
    p(p < 0) = -abs(maxP) * r;        
end


% function Y = compressOutNans(X)
%     nBins = size(X,1);
%     nStims = size(X,3);
%     
%     if nStims == 1
%         [R,L] = meshgrid(1:nBins, 1:nBins);
%         ok_idx = arrayfun(@(r,l) r >= l, R, L);            
% 
%         Y = X(ok_idx);        
%     else
%        
%         [R,L] = meshgrid(1:nBins, 1:nBins);
%         ok_idx = arrayfun(@(r,l) r >= l, R, L);            
%         [ok_r,ok_c] = find(ok_idx);
% 
%         Y = deal( zeros(nnz(ok_idx), nStims) );
%         for n = 1:nStimMax
%             Y(:,n) = X(ok_r, ok_c, n);
%         end        
%     end
% end
% 
% function Y = decompressOutNans(X)
%     N = size(X,1);
%     nBins = (-1 + sqrt(8*N+1))/2;
% 
%     nStims = size(X,2);
%     
%     
%     if nStims == 1
%         [R,L] = meshgrid(1:nBins, 1:nBins);
%         ok_idx = arrayfun(@(r,l) r >= l, R, L);            
% 
%         Y = nan(nBins, nBins);
%         Y(ok_idx) = X;
%         
%     else       
%         [R,L] = meshgrid(1:nBins, 1:nBins);
%         ok_idx = arrayfun(@(r,l) r >= l, R, L);            
%         [ok_r,ok_c] = find(ok_idx);
% 
%         Y = deal( zeros(nnz(ok_idx), nStims) );
%         for n = 1:nStimMax
%             Y(ok_r, ok_c,n) = X(:, n);
%         end
%         
%     end
% end
% 


 
 %     elseif if NStim > 0                
%         error('tmp');
%         nStimMax = 40;
% %         [l_bin, r_bin] = getBestPsthWindow(data0, W0);
% %         new_order = ord(mean(stimPSTH_vals(l_bin:r_bin,:), 1), 'descend'); 
% 
% %         stimPSTH_vals_ordered = stimPSTH_vals(:,new_order);
%         stimPSTH_vals_cum = cummean(stimPSTH_vals_ordered(:,1:nStimMax), 2);
%         
%         [pvals, tstats, areas]  = deal( nan(nBins, nBins, nStimMax) );
%         
%         for nstim = 1:nStimMax
%             cur_psth = stimPSTH_vals_cum(:,nstim);
%             
%             for l_bin = binStart : binEnd-1
%                 for r_bin = l_bin + (minWindowWidth-1) : binEnd;                    
%                     [pvals(r_bin, l_bin, nstim), tstats(r_bin, l_bin, nstim)] = ...
%                         getOspStatVsPsthBinning(stimPSTH_vals_oe, l_bin, r_bin, [nOri, nSp, nPh], cur_psth, statTypes);
%                     areas(r_bin, l_bin, nstim) = getPsthFracArea(cur_psth, l_bin, r_bin);
%                 end
%             end
%         end
%         
%     end





% 
% s = fieldnames(allPsthWindowData_S2);
% cellFieldName_fcn_old = @(Gid, cellId) sprintf('Gid_%04d_cell_%d_0', Gid, cellId);
% cellFieldName_fcn_new = @(Gid, cellId) sprintf('Gid_%04d_cell_%d_st', Gid, cellId);
% for i = 1:length(s)
%     si = s{i};
%     A = sscanf(si, 'Gid_%d_cell_%d');
%     [gid, cellid] = dealV(A);
%     old_name = cellFieldName_fcn_old(gid, cellid);
%     new_name = cellFieldName_fcn_new(gid, cellid);
%     assert( strcmp(old_name, si) );
%     
%     allPsthWindowData_S.(new_name) = allPsthWindowData_S2.(old_name);
% end
% 
% statTypes = {'rep_p', 'rep_tstat',  'rho', 'rho_p', 'r_entropy'};
% 
% structfun(@(s) rmfield(s, {'rho_t', 'r_var'}), allPsthWindowData_S, 'un', 0)
% 


% psthWindowDataFile = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData_st_compressed'];
% psthWindowDataFile2 = [CatV1Path 'MatLabDB_avi' filesep 'psthWindowData_st'];
% 
% S1 = load(psthWindowDataFile);
% cell_names = fieldnames(S1);
% cellFieldName_fcn = @(s) sprintf('Gid_%04d_cell_%d_st', Gid, cellId);            
% 
% for i = 1:length(cell_names)
%     v = S1.(cell_names{i});
%     stat_names = fieldnames(v);
%     
%     ispval = find( cellfun(@(s) ~isempty(strfind(s, '_p')), stat_names ) );
%     for j = ispval
%         st = v.(stat_names{j});
%         st(1:end-2) = pval2NegLogPval(  st(1:end-2)  );
%         
%         v.(stat_names{j}) = st;
%         if any(v.(stat_names{j})(:) == 0)
%             3;
%         end
%     end    
%     v_decomp = structfun(@decompressOutNans, v, 'un', 0);
%     
%     v_sparse = structfun(@(x) sparse(double(x)), v_decomp, 'un', 0);
%     
%     S2.(cell_names{i}) = v_sparse;    
% end
% 
% 

function idxs = idxsOfBinA(A, B)    
    iB = 1; nB = length(B);
    iA = 1; nA = length(A);
    idxs = zeros(1, max(nA, nB));
    ii = 1;
    while (iA <= nA) && (iB <= nB) 
        if strcmp(A{iA}, B{iB})
            idxs(ii) = iA;
            ii = ii+1;
            iA = iA+1;
            iB = iB+1;
        else
            iA = iA+1;
        end        
    end
    idxs(ii:end) = [];        
    if length(idxs) < length(B)
        error('not all elements of B found')
    end
end
        
        
        





