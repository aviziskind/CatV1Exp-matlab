function [psthBins, psthVals, psthStats] = getCellPSTHvals(Gid, cellId)

    persistent allCellPSTHs saveCount matchDB_here
        
    
    opts.nStimMax = 30;
    opts.stim_ordering = 'mean';
    opts.sortingWindow_ms = [20, 150];    
    opts.psthWindow = [-300, 200];
    
    redo_all = 0;
    redo_current = 0;
    saveCountSpacing = 50;
    
    fet_str = iff(curMatchDB, '_DB', ['_' curSortingFeatures('')]);
    cellPSTHs_file = [CatV1Path 'MatLabDB_avi' filesep 'allCellPSTHs' fet_str '.mat'];        
    
    if strcmp(Gid, 'save') && ~isempty(allCellPSTHs)
        save(cellPSTHs_file, 'allCellPSTHs', '-v6');                
        saveCount = 0;
        return;
    end        
        
    if isempty(allCellPSTHs) || (curMatchDB ~= matchDB_here)
        if exist(cellPSTHs_file, 'file') && ~redo_all
            S_file = load(cellPSTHs_file);
            allCellPSTHs = S_file.allCellPSTHs;
        else
            allCellPSTHs = struct;
        end        
        saveCount = 0;
        matchDB_here = curMatchDB;
    end
        
    
    cell_fld_name = sprintf('PSTH_Gid_%d_cellId_%d', Gid, cellId);
    cell_fld_name = strrep(cell_fld_name, '-1', 'n1');

%     if (~isfield(allCellPSTHs, cell_fld_name) || redo_current)
%         allPSTHs = calcBackgroundSpikes(Gid, cellId);
    differentOpts = 0;
    if isfield(allCellPSTHs, cell_fld_name) 
        opts_used = allCellPSTHs.(cell_fld_name).opts;        
        differentOpts = ~isequal(opts, opts_used);
    end

    if (~isfield(allCellPSTHs, cell_fld_name) || redo_current || differentOpts)        
        
        [PSTH_bins, PSTH_vals, PSTH_stats] = calcPSTHforFlashGratingCells(Gid, cellId, opts.psthWindow, opts);        
        PSTH_data.bins = PSTH_bins;
        PSTH_data.vals = PSTH_vals;
        PSTH_data.stats = PSTH_stats;
        PSTH_data.opts = opts;
                        
        allCellPSTHs.(cell_fld_name) = PSTH_data;
        saveCount = saveCount + 1;
        
        if saveCount > saveCountSpacing
            save(cellPSTHs_file, 'allCellPSTHs', '-v6');        
            saveCount = 0;
        end                
    end
        
    S = allCellPSTHs.(cell_fld_name);
    [psthBins, psthVals, psthStats] = deal(S.bins, S.vals, S.stats);   

end

