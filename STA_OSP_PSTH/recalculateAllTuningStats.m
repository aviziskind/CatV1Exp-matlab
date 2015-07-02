function recalculateAllTuningStats(stimType, onlyOnesInDF_flag, onlyCells_flag)

    global GC
    global doOddEvenTrials
    changeInOSPDataFile = 1;
%     global showWorking
%     showWorkingTmp = showWorking;
%     showWorking = false; %#ok<NASGU>
    
%     gratingType = 'flashed';
    
%     [gratingTypeId, gratingType] = curGratingType;
%     if ~exist('gratingType', 'var') || ~ischar(gratingType)
%         gratingType = input('Grating Type? ', 's');
%     else
%     fprintf(['Recalculating stats for ' gratingType ' gratings ...\n'])
%     end
                    
%     if ~exist('allOSPs', 'var')
%         load 
%             
%     s = who('celldata__*');

    if nargin < 1
        fprintf('Recalculating ALL tuning stats (orientation batches, spatial frequencies, flashed gratings...\n');        

%         recalculateAllTuningStats('s', 1, 1);
%         recalculateAllTuningStats('f', 1, 1);
%         recalculateAllTuningStats('o', 1, 1);


%         recalculateAllTuningStats('s', 1, 1);
        recalculateAllTuningStats('f', 1);
        recalculateAllTuningStats('s', 1);
%         recalculateAllTuningStats('o', 1);
% 
%         recalculateAllTuningStats('s', 0, 1);
%         recalculateAllTuningStats('f', 0, 1);
%         recalculateAllTuningStats('o', 0, 1);



        recalculateAllTuningStats('s');

        recalculateAllTuningStats('f');
        recalculateAllTuningStats('o');

        return;
    end

    if ischar(stimType)
        switch stimType
            case 'o',  curGratingType('d'); filename = getFileName('indiv', 'grating_dOr');  stimType_str = 'Orientation Tuning Batches';
            case 's',  curGratingType('d'); filename = getFileName('indiv', 'grating_dSf');  stimType_str = 'Spatial Frequency Tuning Batches';
            case 'f',  curGratingType('f'); filename = getFileName('indiv', 'movie_fg');     stimType_str = 'Flashed Grating Batches';
            otherwise, filename = stimType;
        end
        S_file = load(filename);        
        nFieldNames_orig = length(fieldnames(S_file));
    end
    
    onlyOnesInDF = exist('onlyOnesInDF_flag', 'var') && isequal(onlyOnesInDF_flag, 1);
    onlyCells = exist('onlyCells_flag', 'var') && isequal(onlyCells_flag, 1);
%     onlyCells
    
    if onlyOnesInDF      
        osp_filename = getFileName('osps', '_all', [], struct('windowOffset', 0));
        S_cells = load(osp_filename);
                
        switch stimType
            case 'o',  idx_stim = find(strncmp({S_cells.allCells.stimType}, 'Grating:Orientation', 15)); 
            case 's',  idx_stim = find(strncmp({S_cells.allCells.stimType}, 'Grating:Spatial Freq', 15)); 
            case 'f',  idx_stim = find(strncmp({S_cells.allCells.stimType}, 'Movie:Flashed', 10)); 
        end
        Gids = [S_cells.allCells(idx_stim).Gid];
        cellIds = [S_cells.allCells(idx_stim).cellId];
        subset_str = 'used';
    else               
        allVarnames = fieldnames(S_file);
        [Gids, cellIds] = cellfun(@getGidCellIdFromVarname, allVarnames);
        subset_str = 'ALL';
    end
        
    if onlyCells
        %%
        idx_cells = cellIds > 0 & cellIds < 100;
        Gids = Gids(idx_cells);
        cellIds = cellIds(idx_cells);  
        
%         Gids_plot = [5144];
%         cellIds_plot = [4, 1, 7, 2, 3];
% %         
%      Gids_plot = [701 714 714 716 830 1140 2034 2629 2739 2863 3045 3075 3089 4001 4001 4087 4504 4540 4624 5000 5050 5148 5148 5166 5172 5172 5272]; 
%      cellIds_plot = [ 2 3 6 1 2 1 2 1 3 2 2 8 6 1 3 2 3 4 1 1 1 1 6 2 1 3 2];
%  
% %  
%         idx_use = [];
%         for i = 1:length(Gids)
%             for j = 1:length(Gids_plot)
% %                 if Gids(i) == Gids_plot(j) && cellIds(i) == cellIds_plot(j)
%                 if Gids(i) == Gids_plot(j)
%                     idx_use = [idx_use, i];
%                 end
%             end
%         end
% %         3;
% % %         idx_use = cellfun(@(gc) @([
% % % %         idx = binarySearch(Gids, Gids_plot, 1, 0);
% % %         idx_use = find( any(bsxfun(@eq, Gids(:)', Gids_plot(:)), 1) & any(bsxfun(@eq, cellIds(:)', cellIds_plot(:)), 1)  );
%         Gids = Gids( idx_use );
%         cellIds = cellIds( idx_use );  

    end
    
    
%     S2 = load('flashedGratingCells_GLFcuw8_degree_all.mat');
%     Gids = [S2.allCells.Gid];
%     cellIds = [S2.allCells.cellId];
    
%     cell_gTypeId = flashedOrDrifting(Gids);    
    
%     idx = find(cell_gTypeId == gratingTypeId);
        
    nCells = length(Gids);    
    fprintf('Recalculating %s degree-of-tuning stats for %s (%d cells) \n', subset_str, stimType_str, nCells);
    
%     nRecalced = 0;
    showIndivProgress = 1;

%     if ~showIndivProgress
%         progressBar('init-', nCells, 60);
%     end
    for i = 1:nCells
        
        Gid = Gids(i);
        cellId = double(cellIds(i));
        3;
        sd = siteDataFor('Gid', Gid, 1);
%         varname = getName('celldata', Gid, cellId);
        
        if ~any(sd.cellIds == cellId) && (cellId < 100) && isfield(S_file, varname)
            keyboard;
%             S_file = rmfield(S_file, varname);            
            continue;
        end
           %%     
        GC = [Gid, cellId];        
%         v = S_file.(varname);
        
        if Gid == 1753
            3;
        end

        if changeInOSPDataFile
%             PSTHdata = v.PSTH;
            isflashed = flashedOrDrifting(Gid) == 1;
            if isflashed
                PSTHdata = getPSTHforCell(Gid, cellId);
                LR_bins = PSTHdata.timeWindow_bins;
                windowProfile = PSTHdata.windowProfile;
            else
                PSTHdata = [];
                LR_bins = [1 1];
                windowProfile = [];
            end


%             isflashed = ~isempty(PSTHdata); 
%             if isflashed
%                 LR_bins = v.PSTH.timeWindow_bins;
%                 windowProfile = v.PSTH.windowProfile;
%                 
%                 extraBinOffsets = {}; %{[-1, 1], [1, -1], [-2, 2], [2, -2]};
%                 LR_bins_stimw = v.PSTH.timeWindow_stimw;
%                 windowProfile_stimw = v.PSTH.windowProfile_stimw;
%                 
%             else                
%                 LR_bins = [1 1];
%                 windowProfile = [];
%                 extraBinOffsets = [];                
%             end
%             redo_flag = 735878.574639; %735830.154287; % 735819.711500; %735384.066367; %; %sprintf('%.6f', now)
%             redo_flag = 735830.154287; % 735819.711500; %735384.066367; %; %sprintf('%.6f', now)

            redo_flag = 735993.628037;
            
            if showIndivProgress
                tic;
                fprintf('(%d/%d) Gid = %d. cellId = %d. ', i, nCells, Gid, cellId);
            end
            doOddEvenTrials = false;
            if ~isempty(PSTHdata)
                fprintf('[%d-%d]', LR_bins);
            end
            tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'}, redo_flag);
            
%             if ~isfield(tuningStats.oriStats_ss.error_jack, 'F1oDC')
%                 tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'}, 1);
%             end
            
%             if (length(fieldnames(tuningStats)) == 4)
%                 tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'}, 1);
%             end
            
            if isflashed && 0
                doOddEvenTrials = false;
%                 fprintf('[%d-%d]', LR_bins_stimw);
%                 tuningStats = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins_stimw(1), LR_bins_stimw(2), windowProfile_stimw, {'tuningStats'}, redo_flag);
                if between(cellId, 0, 100)
                    for j = 1:length(extraBinOffsets)
                        
                        LR_bins_j = LR_bins + extraBinOffsets{j};
                        if diff(LR_bins_j) < 2
                            continue
                        end
                        windowProfile_j = v.PSTH.vals(LR_bins_j(1):LR_bins_j(2));
                        fprintf('[%d-%d]', LR_bins_j);
                        tuningStats_j = getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins_j(1), LR_bins_j(2), windowProfile_j, {'tuningStats'}, redo_flag);
                    end
                end
            end
            
            if showIndivProgress
                toc;
            end
        else                    
            tuningStats = calcDegreeOfTuningStats( v.OSP.R_full, v.OSP.bckgSamples, v.Gid);        
        end
%         v.OSP.stats.tuningStats = tuningStats;
%         S_file.(varname) = v;        
        
%         progressBar;
    end
%     if ~showIndivProgress
%         progressBar('done');
%     end    
    
    S_file = orderfields(S_file);
    fprintf('Recalculated stats for %d cells ... \n', nCells);
        
    nFieldNames_now = length(fieldnames(S_file));
    assert(nFieldNames_now >= nFieldNames_orig);
    
    save(filename, '-struct', 'S_file', '-v6');
    
    if changeInOSPDataFile
        %%
        dbGetCellSpkStimHists('save');  
        getOspDataForPsthWindow('save');
    end
    
end
    
% end


%{
for i = 1:length(idx), 
    ss{i} = calculatePSTH_STAs_OSP_ForOneCell(S.allCells(idx(i)).Gid, S.allCells(idx(i)).cellId);
end

for i = 1:length(idx), 
    ss{i} = calculatePSTH_STAs_OSP_ForOneCell(S.allCells(idx(i)).Gid, S.allCells(idx(i)).cellId);
    figure(i);
    imageOSP( s.OSP.R, 'mean:ph');
    title(sprintf('Group %d, cell %d. (bckg = %.1f Hz)', s.Gid, s.cellId, mean(decompress(s.OSP.bckgSamples))))
end
%}