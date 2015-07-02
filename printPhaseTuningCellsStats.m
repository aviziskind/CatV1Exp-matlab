function printPhaseTuningCellsStats

    printCellSelectionCriteriaStats = 1;
    printStatsOfGoodCells = 0;
    
    % before doing this: run generateGratingCellsDatafile with 'all' selected

    curCmpType('phase');        
    gratingType = curGratingType('');  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;    
    timeWindow_str = curTimeWindow('');
%     [cmpType, cmpType_s] = curCmpType;        
        

    
    fprintf(' ***************************** \n ****** %s GRATINGS (timeWindow: %s ****** \n *****************************', upper(curGratingType('')));

    
    if printCellSelectionCriteriaStats
        
        min_ID = 8;
        
        alpha = .001;
        log_alpha = 3;
        minJackCC = -1; %0.5;
        minRsqr_oe = 0.5;
        minRsqr_fit = 0.4;
        
        % Cell selection criteria statistics
        
        if strcmp(gratingType, 'flashed')
            allTimeWindows = {'best', [30 60], [60 90]};
        elseif strcmp(gratingType, 'drifting')
            allTimeWindows = {'best'};
        end
        
        frac_pct = @(i,n) sprintf('%3d/%3d (%.1f%%)', i, n, i/n*100);
        frac_pct_tf = @(tf) frac_pct( nnz(tf), length(tf) );
        allMeasureNames = {};
        for ti = 1:length(allTimeWindows)
            
             curTimeWindow(allTimeWindows{ti})
            fprintf(' ***************************** \n  (timeWindow: %s ****** \n *****************************', curTimeWindow('')); 
    
           
            ospDatafile = getFileName('osps', '_all');  % use file which includes non-reproducible/selective cells
            % (but still has narrowed down to 1 group per site)
            S1 = load(ospDatafile);
            allCells = S1.allCells;
            clear S1;
            %     nUnits = length(allCells);
            
            idx_cells = find([allCells.cellId] > 0);
            allCells = allCells(idx_cells);
            nCellsTot = length(idx_cells);
            
            switch gratingType,
                case 'flashed'
                    
                case 'drifting',
                    allCells = allCells( arrayfun(@(s) strncmp(s.stimType, 'Grating:Spatial Freq', 19), allCells ) );
            end
            
            fprintf('\n\n Cell Profile Selection criteria\n');
            IDs = nestedFields(allCells, 'spkFeatures', 'IsolationDistance');            
            idx_wellIsolated = IDs > min_ID;
            nCellsIsolated = nnz( idx_wellIsolated );
            fprintf([' Phase : Well-isolated : %s \n'], frac_pct(  nCellsIsolated , nCellsTot) );
            
            
            
            allCells = allCells(idx_wellIsolated);
            nCells = length(allCells);
            %     oriUsableCells = o
%             stats = [allCells.stats];
%             windowStats = [allCells.windowStats];

            windowStats = nestedFields(allCells, 'windowStats', -1);
           
            
            
            idx_phase_reproducible = [windowStats.cc_p] > log_alpha;
            idx_selected = [allCells.cell_ok];

            
    %         goodCells = allCells(idx_selected);
    %         nGoodCells = length(goodCells);        
            nCells_used = nnz(idx_phase_reproducible);

%             f1odc_field = iff(strcmp(gratingType, 'flashed'), 'F1oDC_maxR_avP', 'F1oDC_maxR_avP_sm');
            idx_simpleCells = [allCells.F1oDC_maxR_avP_sm] > 1;
            idx_complexCells = ~idx_simpleCells;

        %     ori_nAccepted = nnz([oriCellStats_si    

            fprintf('\n\n Cell Profile Selection criteria\n');
            fprintf(['  Phase : Reproducible : %s \n'], frac_pct_tf(  idx_phase_reproducible ) );
            fprintf([' Phase : Total currently used : %s \n']', frac_pct_tf(  idx_selected  ) );

            blankMIDstruct = struct('MID', [], 'MID_select', [], 'MID_fit', [], 'jackCC', nan, 'rsqr_fit', nan, 'p', [], 'rsqr_oe', nan, 'pval_oe', nan);

            if strcmp(gratingType, 'flashed');
                %%
                idx_reproducible_andHaveMID = idx_selected; % & arrayfun(@(s) isfield(s.MIDdata, 'MID'), allCells)';

                MID_flds = {'MID', 'MID_odd', 'MID_even'};
                for i = 1:length(allCells)
                    haveMIDs_i = 1;
                    for j = 1:length(MID_flds)
                        if ~isfield(allCells(i).MIDdata, MID_flds{j})                        
                            allCells(i).MIDdata.(MID_flds{j}) = blankMIDstruct;
                            haveMIDs_i = 0;
                        end
                    end
                    idx_reproducible_andHaveMID(i) = idx_selected(i) & haveMIDs_i;
    %                 if ~isfield(allCells(i).MIDdata, 'rsqr_oe')
    %                     allCells(i).MIDdata.rsqr_oe = nan;
    %                     allCells(i).MIDdata.pval_oe = nan;
    %                 end
                end

%                 goodCellsWithMID = allCells(idx_reproducible_andHaveMID);            


                
    %             MID_rep_use = 'jackCC';
                MID_rep_use = 'rsqr_oe';
                minRsqr_fit = 0.5;

                all_MIDdata = nestedFields(allCells, 'MIDdata', 1);
%                 nCellsWithMID = length(all_MIDdata);
    %             jackCC = [all_MIDdata.jackCC];
    %             rsqr   = [all_MIDdata.rsqr];
%                 assert(isequal (idx_used_simpleCells, [goodCellsWithMID.F1oDC_maxR_avP_sm] > 1) );
                nSimpleCells_used = nnz( idx_selected & idx_simpleCells);
                nComplexCells_used = nnz( idx_selected & ~idx_simpleCells );
    %             idx_MID_rep = [all_MIDdata.jackCC].^2 > minJackCC;

                idx_simpleCells_withMID = idx_simpleCells(idx_reproducible_andHaveMID);
                idx_complexCells_withMID = ~idx_simpleCells(idx_reproducible_andHaveMID);
                
                nSimpleCellsWithMID = nnz(idx_simpleCells_withMID);
                nComplexCellsWithMID = nnz(idx_complexCells_withMID);
                
                all_rsqr_oe = [all_MIDdata.rsqr_oe];            

                jackCC_odd = nestedFields(all_MIDdata, 'MID_odd', 'jackCC', 1).^2;
                jackCC_even = nestedFields(all_MIDdata, 'MID_even', 'jackCC', 1).^2;            
                all_smaller_jackCC_fit = min( jackCC_odd, jackCC_even);

                if strcmp(MID_rep_use, 'jackCC')
                    idx_MID_rep = all_smaller_jackCC_fit > minJackCC;    
                    min_r = minJackCC;
                elseif strcmp(MID_rep_use, 'rsqr_oe')
                    idx_MID_rep = all_rsqr_oe > minRsqr_oe;            
                    min_r = minRsqr_oe;
                end

%%
                rsqr_fit_odd = nestedFields(all_MIDdata, 'MID_odd', 'rsqr_fit', 1);
                rsqr_fit_even = nestedFields(all_MIDdata, 'MID_even', 'rsqr_fit', 1);            
                all_smaller_rsqr_fit = min( rsqr_fit_odd,  rsqr_fit_even);
                idx_MID_good_fit = all_smaller_rsqr_fit > minRsqr_fit;

                fprintf('\n\n ALL MID Selection criteria\n');
                fprintf([' MID : Reproducible (%s > %.2f): %s \n'], MID_rep_use, min_r, frac_pct( nnz(idx_selected & idx_MID_rep) , nCells_used) );            
                fprintf([' MID : Reproducible & Good Fit (%s > %.2f, rsqr > %.2f): %s  [ %s ] \n'], MID_rep_use, min_r, minRsqr_fit, ...
                    frac_pct( nnz( idx_MID_rep & idx_MID_good_fit), nCells_used ), ...
                    frac_pct( nnz( idx_MID_rep & idx_MID_good_fit), nCells_used ) ...
                );

                fprintf('\n\n Simple cell MID Selection criteria\n');
                fprintf([' MID : # Simple cells (F1/DC > 1): %s \n'],   frac_pct(  nSimpleCells_used, nCells_used) );            
                fprintf([' MID : Reproducible (%s > %.2f): %s \n'], MID_rep_use, min_r, frac_pct(  nnz( idx_selected & idx_simpleCells & idx_MID_rep ),  nnz( idx_selected & idx_simpleCells )  ) );            
                fprintf([' MID : Reproducible & Good Fit (%s > %.2f, rsqr > %.2f): %s  [%s]\n'], MID_rep_use, min_r, minRsqr_fit, ...
                    frac_pct( nnz(  idx_selected & idx_simpleCells & idx_MID_rep & idx_MID_good_fit), nSimpleCells_used),...
                    frac_pct( nnz(  idx_selected & idx_simpleCells & idx_MID_rep & idx_MID_good_fit), nnz( idx_simpleCells & idx_MID_rep ) ));
                3;

                fprintf('\n Complex cell MID Selection criteria\n');
                fprintf([' MID : # Complex cells (F1/DC > 1): %s \n'],   frac_pct(  nComplexCells_used, nCells_used) );            
                fprintf([' MID : Reproducible (%s > %.2f): %s \n'], MID_rep_use, min_r, frac_pct(  nnz( idx_selected & idx_complexCells & idx_MID_rep ),  nnz( idx_selected & idx_complexCells ) ) );            
                fprintf([' MID : Reproducible & Good Fit (%s > %.2f, rsqr > %.2f): %s  [%s]\n'], MID_rep_use, min_r, minRsqr_fit, ...
                    frac_pct(  nnz( idx_complexCells & idx_MID_rep & idx_MID_good_fit), nComplexCells_used),...
                    frac_pct(  nnz( idx_complexCells & idx_MID_rep & idx_MID_good_fit), nnz( idx_complexCells & idx_MID_rep ) ));

                3;



            else

                [idx_MID_rep, idx_MID_good_fit] = deal(false);
            end

            %%
            stats_name = sprintf('stats%s', curTimeWindow(''));
            allStats.(stats_name) = struct('n_total_isol', nCellsTot, 'n_total', nCells, ... 'n_reproducible', nnz( idx_phase_reproducible ), ...
                'n_used', nnz( idx_selected ), ...
                'n_used_simple', nnz( idx_selected & idx_simpleCells ), ...
                ...
                'n_MID_rep', nnz( idx_selected & idx_MID_rep ), ...
                'n_MID_rep_fit', nnz( idx_selected & idx_MID_rep & idx_MID_good_fit ), ...
                ...
                'n_MID_rep_simple', nnz( idx_selected & idx_simpleCells & idx_MID_rep ), ...
                'n_MID_rep_fit_simple', nnz( idx_selected & idx_simpleCells & idx_MID_rep & idx_MID_good_fit ) ...
                ...
                ... 'pct_used_spf', nnz( idx_spfSelected ) / spf_nCells * 100 ...            
                 ...
            );        

    %         data_crit_S.S_criteria_spf = S_criteria_spf;
            %%
            allMeasureNames = [allMeasureNames, stats_name]; %#ok<AGROW>
            data_crit_S.(stats_name) = allStats.(stats_name);
            
%             data_crit_S.timeWindow = curTimeWindow('');
        end
        data_crit_S.allMeasureNames = allMeasureNames;
        allColumns = fieldnames(allStats.(stats_name));
        data_crit_S.columns = allColumns; %fieldnames(data_S.(allMeasureNames{i}));
        data_crit_S.cmpType = curCmpType('');
        data_crit_S.gratingType = curGratingType('');
        
        criteria_fn = getFileName('cellSelectedStats', [], [], struct('gratingType', gratingType, 'cmpType', 'phase', 'timeWindow', 'best'));
        save(criteria_fn, '-struct', 'data_crit_S');    
        
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%     windowStats = [stats.allWindowStats];
            
            
    
    if printStatsOfGoodCells
        
        curTimeWindow('best');
        ospDatafile = getFileName('osps', '_all');  % use file which includes non-reproducible/selective cells
        % (but still has narrowed down to 1 group per site)
        S1 = load(ospDatafile);
        allCells = S1.allCells;
        clear S1;
        %     nUnits = length(allCells);

        idx_cells = find([allCells.cellId] > 0);
        allCells = allCells(idx_cells);

        switch gratingType,
            case 'flashed'

            case 'drifting',
                allCells = allCells( arrayfun(@(s) strncmp(s.stimType, 'Grating:Spatial Freq', 19), allCells ) );
        end
        %     oriUsableCells = o
        windowStats = nestedFields(allCells, 'windowStats', -1);
        
        
        
        fprintf('\n\n  Parameter      |    Mean     |      Std    |    Median   |     P25     |     P75      | N\n');

    %     statNames = {'w_ori_global_cell', 'w_ori_local_cell', 'w_ori_global_MU', 'w_ori_local_MU', 'dsi_cell', 'dsi_MU', 'spf_width', 'spf_pref'};
        
        ptcStatNames = {'cc', 'dphi', 'cc-MID', 'cc-MID-fit'};     
        MIDStatNames = {'cc-MID', 'cc-MID-fit'};
        
        allStatNames = [ptcStatNames, MIDStatNames];    
        
        all_cells_ok = [oriCellStats_si.cellOK];
        all_cellStats_ok = oriCellStats_si(ori_cells_ok);                    

        
        for i = 1:length(allStatNames)
            statName = allStatNames{i};

            if strcmp(gratingType, 'drifting') && any(strcmp(statName, MIDStatNames))
                continue;
            end
            
            switch statName
                case ptcStatNames{1},  
                    fprintf('*** PHASE TUNING CURVE  *** \n');
                    
                case MIDStatNames{1},  
                    fprintf('*** MID''s *** \n');                    
            end
            
            nSpont = 1;
                        
            [vals_mean, vals_std, vals_median, vals_P25, vals_P75, N] = deal( cell(1, nSpont) );
            for j = 1:nSpont
                vals = [all_cellStats_ok.(statName)];
                [vals_mean{j}, vals_std{j}, vals_median{j}, vals_P25{j}, vals_P75{j}, N{j}] = getMeanStdStats(vals);  
            end        


            w = num2str( switchh( nSpont, [1 2], [11, 5]));
            if strncmp(statName, 'dphi', 3) 
                fcode = ['%' w '.1f'];
            else % strncmp(statName, 'DSI', 3) || any(strcmp(statName, suppStatNames)) || any(strcmp(statName, spfStatNames))
                fcode = ['%' w '.2f'];
            end

            if nSpont == 2
                fcode_tmp = ['| ' fcode '/' fcode ' '];
                n_tmp = '%d/%d';
            else
                fcode_tmp = ['| ' fcode  ' '];
                n_tmp = '%d';
            end

            str_template = ['%16s ' repmat( fcode_tmp, 1, 5) ' | ' n_tmp ' \n'];        
            fprintf(str_template, ...
                allStatNames{i}, vals_mean{:}, vals_std{:}, vals_median{:}, vals_P25{:}, vals_P75{:}, N{:});
                
            
        end
    end

    

end







function [x_mean, x_std, x_median, x_p25, x_p75, n] = getMeanStdStats(x)
     x_mean = nanmean(x);   
     x_std = nanstd(x);
     x_median = nanmedian(x);
     x_p25 = prctile(x,25);
     x_p75 = prctile(x,75);
     n = nnz(~isnan(x));    
end
    
    


% %%
% Sd = load('driftingGratingCells_GLFcuw8_degree_all.mat');
% idx_ori_cells_d = find(arrayfun(@(s) s.cellId > 0 & length(s.ori)>=30, Sd.allCells));
% Sd_ori_stats = nestedFields(Sd.allCells(idx_ori_cells_d), 'stats', 'tuningStats', 'oriStats_si');
% idx_ok_d = find([Sd_ori_stats.cellOK]);
% Sd_ori_stats_params = nestedFields(Sd_ori_stats(idx_ok_d), 'oriParams');
% 
% [mean([Sd_ori_stats_params.B]), std([Sd_ori_stats_params.B])]
% %%
% Sf = load('flashedGratingCells_GLFcuw8_degree_all.mat');
% idx_ori_cells_f = find(arrayfun(@(s) s.cellId > 0 & length(s.ori)>=30, Sf.allCells));
% Sf_ori_stats = nestedFields(Sf.allCells(idx_ori_cells_f), 'stats', 'tuningStats', 'oriStats_si');
% idx_ok_f = find([Sf_ori_stats.cellOK]);
% Sf_ori_stats_params = nestedFields(Sf_ori_stats(idx_ok_f), 'oriParams');
% 
% [mean([Sf_ori_stats_params.B]), std([Sf_ori_stats_params.B])]