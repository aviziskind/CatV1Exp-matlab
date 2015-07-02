function generateExcelSpreadsheets

%     if nargin < 1
% %         whichData = 'cellStats';
%         whichData = 'pairStats';
%     end    
	opt.separateCellsForFlashedAndDrifting = 1;
    opt.addExtraRowForDF = 1;
    opt.includeComments = 0;
    doDegree = 1;
        doRegular_ss = 0;
        doRegular_si = 0;

        doSamePen_ss = 0;
        doSamePen_si = 0;

        doSameAnimal_ss = 0;
        
        doBPWA_ss = 0;

        doSimpComplex_ss = 1;
        doSimpComplex_si = 0;

        doAlignedAntialigned = 0;
        
        doOddEvenTests = 1;
            doOE_same = 1;
            doOE_diff = 1;
            doOE_same_simpComplex = 1;
            doOE_diff_simpComplex = 1;
    
    doPhase = 0;
%     cmpType = 'degree';
%     cmpType = 'phase';
%     curPreserveSimpleComplex(0);
%     curSamePen(0);
    redoAll = 1;    
%         redoBefore = 735965.523747; % sprintf('%.6f', now)
        redoBefore = 735982.533527; % sprintf('%.6f', now)

        
    if redoAll
        redoBefore = inf;
    end
%     doCellStats = 1;
%     doPairStats = 1;

%                                  subtSpont  preserveSC   bccType
%     allDegreePairStatsCombos = {      0          0          0     ;  ... standard degree of tuning test.
%                                       0          0          1     ;  ... only same penetration (test for spf)
%                                       0          1          0     ;  ... for simple/complex clustering analysis
%                                       1          0          0     };  ... for spontaneous subtracted measures.
%     allPhasePairStatsCombos = {       0          0          0     }; ... phase tuning. (no between site, and no Wrcc, so preserveSC and bccType are irrelevant)                   

    
    if doDegree            
        cmpType = 'degree';
        timeWindow = 'best';
        
        
        3;
        
        if doRegular_ss
             % spont-subtracted measure
            subtractSpont = 1; bccType = 'full'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'aa';
%             generateSpreadsheet('cellSelectedStats',       cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
%             generateSpreadsheet('cellStats',                cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
%             generateSpreadsheet('cellDistribs',                cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
%             generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
%             generateSpreadsheet('pairStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);            

            generateSpreadsheet('scCellStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt); % use preserveSC for Simple/Complex ;
            generateSpreadsheet('scPairStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt); % use preserveSC for Simple/Complex ;
            3;
%             generateSpreadsheet('penDistStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
        end
        
        if doRegular_si
            % standard measure:
            subtractSpont = 0; bccType = 'full'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'aa';
            generateSpreadsheet('cellStats',                cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
            generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);            
        end
        

      
        if doSamePen_ss
            % same penetration only
            subtractSpont = 1; bccType = 'samePen'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'aa';
%             generateSpreadsheet('cellStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore);
            generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
        end     
        
        if doSamePen_si
            % same penetration only
            subtractSpont = 0; bccType = 'samePen'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'aa';
%             generateSpreadsheet('cellStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore);
            generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
        end
        
        if doSameAnimal_ss
            % same penetration only
            subtractSpont = 1; bccType = 'sameAnimal'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'aa';
            %             generateSpreadsheet('cellStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore);
            generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
        end
        
        if doBPWA_ss
            % same penetration only
            subtractSpont = 1; bccType = 'diffPen-sameAnimal'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'aa';
            %             generateSpreadsheet('cellStats',              cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore);
            generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
        end
        
       if doSimpComplex_ss
            % simple/complex clustering analysis
            subtractSpont = 1; bccType = 'full'; preserveSC = 1; preserveAB = 0;  degree_oe_mode = 'aa';
            generateSpreadsheet('scCellStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
            generateSpreadsheet('scPairStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
            3;
       end     
        
        if doSimpComplex_si

            % simple/complex clustering analysis
            subtractSpont = 0; bccType = 'full'; preserveSC = 1; preserveAB = 0;  degree_oe_mode = 'aa';
            generateSpreadsheet('scCellStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
            generateSpreadsheet('scPairStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
        end       

  
        
        if doAlignedAntialigned
            subtractSpont = 1; bccType = 'full'; preserveSC = 0; preserveAB = 1;  degree_oe_mode = 'aa';
            opt_onlyDG = opt;
            opt_onlyDG.gratingTypes = {'drifting'};
            generateSpreadsheet('pairStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt_onlyDG);                        
        end
        
        
        
        if doOddEvenTests
            
            if doOE_same                
                subtractSpont = 1; bccType = 'full'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'oe_same';
%                 generateSpreadsheet('cellStats',                cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
                generateSpreadsheet({'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
            end
            
            if doOE_diff
                subtractSpont = 1; bccType = 'full'; preserveSC = 0; preserveAB = 0;  degree_oe_mode = 'oe_diff';
%                 generateSpreadsheet('cellStats',                cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);
                generateSpreadsheet({'pairStats'}, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);               
            end
            
            if doOE_same_simpComplex 
                subtractSpont = 1; bccType = 'full'; preserveSC = 1; preserveAB = 0;  degree_oe_mode = 'oe_same';
                generateSpreadsheet('scPairStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);                
            end
            
            if doOE_diff_simpComplex
                subtractSpont = 1; bccType = 'full'; preserveSC = 1; preserveAB = 0;  degree_oe_mode = 'oe_diff';
                generateSpreadsheet('scPairStats', cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, degree_oe_mode, redoBefore, opt);                                
            end

            
        end
        
        
    end
    
        
    if doPhase
        cmpType = 'phase';
        phase_oe_mode = 'oe_diff_k12';
               
%         generateSpreadsheet('cellStats', cmpType, [], [], [], timeWindow, redoBefore);
%         timeWindows = {'best', [29, 62], [58, 91]};
        timeWindows = {'best'};
        for i = 1:length(timeWindows)
            timeWindow = timeWindows{i};
            timeWindow_str = iff(ischar(timeWindow), timeWindow, sprintf('[%d-%d]', timeWindow));
            fprintf('Phase tuning spreadsheets: Doing pair diffs & pair stats for timeWindow : %s\n', timeWindow_str);
%             generateSpreadsheet({'cellSelectedStats'},     cmpType, [], [], [], [], timeWindow, phase_oe_mode, redoBefore, opt);
% %             generateSpreadsheet('pairDiffs',               cmpType, [], [], [], [], timeWindow, redoBefore);
%             generateSpreadsheet('pairStats', cmpType, [], [], [], [], timeWindow, redoBefore, );
            generateSpreadsheet({'pairDiffs', 'pairStats'}, cmpType, [], [], [], [], timeWindow, phase_oe_mode, redoBefore, opt);
        end                
        
    end

    
    
end

function generateSpreadsheet(statType, cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, redoBefore, opt)
    	
	suffix_str = iff(opt.separateCellsForFlashedAndDrifting, '', '_together');
	
    excelFile = [CatV1Path 'tables_' cmpType suffix_str '.xls'];    
        
    redoMatFiles     = 1;% && strcmp(statType, 'cellStats');
    makeSpreadSheets = 1 && 1;
	
    if isfield(opt, 'gratingTypes')
        gratingType_vals = opt.gratingTypes;
    else
        gratingType_vals = {'drifting', 'flashed'}; % = drifting / flashed .
    end
    opt.doDrifting = any(strcmp(gratingType_vals, 'drifting'));
    opt.doFlashed  = any(strcmp(gratingType_vals, 'flashed'));
    
   
    if redoMatFiles        
                
         for g = 1:length(gratingType_vals)
            gratingStr = gratingType_vals{g};
            if ~iscell(statType)
                statType = {statType};
            end                    
            statFilenames = cellfun(@(st) getStatFileName(cmpType, gratingStr, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, st), statType, 'un', 0);
            statFileNames_noPath = cellfun(@extractFilename, statFilenames, 'un', 0);
            statTypes_str = cellstr2csslist(statType);
            statFilenames_str = cellstr2csslist(statFileNames_noPath);    

            switch cmpType
                case 'degree', fprintf('Degree %s : SS=%d. Bcc=%s. pSC=%d. pAB=%d. OE = %s. Grating=%s.  ... ', statTypes_str, subtractSpont, bccType, preserveSC, preserveAB, oe_mode, gratingStr);
                case 'phase',  fprintf('Phase %s : Grating %s ... ', statTypes_str, gratingStr);
            end

            haveFiles = cellfun(@(stat_file) exist(stat_file, 'file') && ~fileOlderThan(stat_file, redoBefore), statFilenames);

            if all(haveFiles)
                fprintf('  -> Already have file(s) : %s \n', statFilenames_str);
            else
                tic; 
                curCmpType(cmpType); curGratingType(gratingStr);
                curSubtractSpont(subtractSpont); curPreserveSimpleComplex(preserveSC);  curBccType(bccType);
                curTimeWindow(timeWindow); curPreserveAligned(preserveAB);
                

                switch cmpType

                    case 'degree', 
                        curDegreeOEmode(oe_mode);
                        switch statTypes_str
                            case 'cellSelectedStats', printDegreeTuningCellsStats(1);
                            case 'cellStats',  printDegreeTuningCellsStats(2);
                            case 'cellDistribs', printDegreeTuningCellsStats([1 2 3]);
                                
                            case 'pairDiffs',  printDegreeOfTuningComparisons(1);
                            case 'pairStats',  printDegreeOfTuningComparisons(2);  
                            case 'pairDiffs, pairStats', printDegreeOfTuningComparisons([1 2]);                            
                            case {'scCellStats', 'scPairStats'}, 
                                               printDegreeOfTuningComparisons(3);
                            case 'penDistStats', printDegreeOfTuningComparisons(4);
                        end
                        

                    case 'phase',                
                        curPhaseOEmode(oe_mode);
                        switch statTypes_str
                            case 'cellSelectedStats', printPhaseTuningCellsStats(1);
                            case 'cellStats', printPhaseTuningCellsStats(1);
                            case 'pairDiffs', printPhaseTuningComparisons(1);
                            case 'pairStats', printPhaseTuningComparisons(2);
                            case 'pairDiffs, pairStats', printPhaseTuningComparisons([1 2]);                            
                        end
                end
                toc;

                madeFiles = cellfun(@(stat_file) exist(stat_file, 'file') && (~fileOlderThan(stat_file, redoBefore) || isinf(redoBefore)), statFilenames);

%                         madeFile = exist(statFilename, 'file') && (~fileOlderThan(statFilename, redoBefore) || isinf(redoBefore));
                assert(all(madeFiles));
                fprintf('  -> Generated files : %s\n', statFilenames_str)
            end
         end
                
                          
        
    end        
    
    if makeSpreadSheets
        if ~iscell(statType)
            statType = {statType};
        end
        for i = 1:length(statType)
            generateThisSpreadSheet(cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, statType{i}, excelFile, opt)
        end
%         switch cmpType        
% 
%             case 'degree',                     
%                 %%            
%                 for subtractSpont_i = subtractSpont_vals                            
%                     generateThisSpreadSheet(cmpType, subtractSpont, bccType, preserveSC, statType, excelFile)
%                 end        
% 
%             case 'phase'            
%                 
%         end
    end
    
end




function generateThisSpreadSheet(cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, statType, excelFile, opt)
    switch statType
        case {'cellStats', 'cellDistribs', 'pairDiffs', 'pairStats'},      preserveSC = 0;
        case {'simpleComplexCells', 'simpleComplexPairs'}, preserveSC = 1;
    end    
    if opt.doDrifting
        [filename_drifting, opt_d] = getStatFileName(cmpType, 'drifting',  subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, statType);
        Sd = load(filename_drifting);        
        confirmCorrectFile(Sd, opt_d);
    else
        Sd = struct;
    end
    
    if opt.doFlashed
        [filename_flashed,  opt_f]  = getStatFileName(cmpType, 'flashed',   subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, statType);    
        Sf = load(filename_flashed);
        confirmCorrectFile(Sf, opt_f);
    else
        Sf = struct;
    end
    
    str_length_max  = 100;

    
    spreadData = {};
    %% Stage 1: if have columns & measures(rows)
    if (isfield(Sf, 'columns') && ~isempty(Sf.columns)) || (isfield(Sd, 'columns') && ~isempty(Sd.columns))
    
        allColumnNames = [{'PropertyName'}, Sf.columns{:}];    
        columns_toRemove_start = {'Wcc', 'Bcc', 'Wrcc', 'Brcc', 'tmp_'};    
        columns_toRemove_any = {'__rand'};    
        
        stringContains = @(s, substrs)  cellfun(@(substr) ~isempty(strfind(s, substr)), substrs); % does the string contain the following substrings (outputs a boolean for each substring)
        
        idx_col_remove = cellfun(@(s)  any(strncmp(s, columns_toRemove_start, 3)) || any( stringContains(s, columns_toRemove_any)), allColumnNames); 
        allColumnNames(idx_col_remove) = [];

        allMeasureNames_d = Sd.allMeasureNames(~cellfun(@isempty, Sd.allMeasureNames));
        allMeasureNames_f = Sf.allMeasureNames(~cellfun(@isempty, Sf.allMeasureNames));
        if length(allMeasureNames_d) > length(allMeasureNames_f)
            allMeasureNames = uniqueInOrder([allMeasureNames_d, allMeasureNames_f]);
        else
            allMeasureNames = uniqueInOrder([allMeasureNames_f, allMeasureNames_d]);
        end
        %%
        % fix order if messed up:
        idx_maxMin = find( strncmp(allMeasureNames, 'maxMinFracR', 9) );
        idx_maxR1xR2 = find( strncmp(allMeasureNames, 'maxR1xR2', 7) );
        if ~isempty(idx_maxMin) && ~isempty(idx_maxR1xR2)
            other_idxs = setdiff(1:length(allMeasureNames), [idx_maxMin, idx_maxR1xR2]);
            allMeasureNames = allMeasureNames([idx_maxMin, other_idxs, idx_maxR1xR2]);
        end
            %%
        measures_toRemove = {'D_aligned_pair', 'D_F1pair'};    
        idx_ms_remove = cellfun(@(s) any(strcmp(s, measures_toRemove)), allMeasureNames); allMeasureNames(idx_ms_remove) = [];

        %%
        nCol_src = length(allColumnNames);
        if opt.separateCellsForFlashedAndDrifting
            if opt.addExtraRowForDF
                colNamesExpanded_name = cellfun(@(cn) {cn, ''}, allColumnNames(2:end), 'un', 0);
                colNamesExpanded_DF   = cellfun(@(cn) {'(d)', '(f)'}, allColumnNames(2:end), 'un', 0);
                allColumnNames_dst = [allColumnNames(1), [colNamesExpanded_name{:}]];
                allColumnNames_DF  = ['Grating',         [colNamesExpanded_DF{:}]];
            else                        
                colNamesExpanded = cellfun(@(cn) {[cn '(d)'], [cn '(f)']}, allColumnNames(2:end), 'un', 0);
                allColumnNames_dst = [allColumnNames(1), [colNamesExpanded{:}]];
            end
            
            nCol_dst = (nCol_src-1)*2+1;
        else
            allColumnNames_dst = allColumnNames;
            nCol_dst = nCol_src;
        end
        allColumnNames_dst = cellfun(@(cn) strrep(cn, 'vals_', ''), allColumnNames_dst, 'un', 0);
        allMeasureNames(strCcmp(allMeasureNames, {'D_F1pair_spf', 'D_F1pair_ori', 'D_alignedPair'})) = [];
        allMeasureNames(strCcmp(allMeasureNames, {'D_aligned_pair_Wcc', 'D_aligned_pair_Bcc',   'D_F1pair_ori_Wcc', 'D_F1pair_ori_Bcc',   'D_F1pair_spf_Wcc', 'D_F1pair_spf_Bcc'})) = [];
        
        nRow = length(allMeasureNames);
        spreadData = cell(nRow+1, nCol_dst);

        spreadData(1,:) = allColumnNames_dst;
        if opt.addExtraRowForDF
            spreadData(2,:) = allColumnNames_DF;
            row_offset = 2;
        else
            row_offset = 1;
        end
        
        for r = 1:nRow
            measureName = allMeasureNames{r};
            
            for c = 1:nCol_src-1
                spreadData{r+row_offset, 1} = measureName;
                column_name = allColumnNames{c+1};

                have_F = isfield(Sf, measureName) && ~excludeField(measureName, 'flashed');
                have_D = isfield(Sd, measureName) && ~excludeField(measureName, 'drifting');
                if have_D                
                    data_d = Sd.(measureName).(column_name);
                    data_d_str = format_data(data_d, measureName, column_name, statType, opt);
                    if length(data_d_str) > str_length_max 
                        error('string too long')
                    end
                end

                if have_F
                    data_f = Sf.(measureName).(column_name);
                    data_f_str = format_data(data_f, measureName, column_name, statType, opt);
                    if length(data_f_str) > str_length_max
                        error('string too long')
                    end
                end

                if opt.separateCellsForFlashedAndDrifting
                    if have_D
                        spreadData{r+row_offset, 2*c} = add_colon(data_d_str);
                    end
                    if have_F
                        spreadData{r+row_offset, 2*c+1} = add_colon(data_f_str);
                    end
                    if ~have_D && ~have_F
                        error('no data?!');
                    end

                else
                    if have_F && have_D
                        data_str = [data_d_str ' / ' data_f_str];
                    elseif have_D
                        data_str = [data_d_str ';'];
                    elseif have_F
                        data_str = [data_f_str ';'];
                    else
                        error('no data?!');
                    end
                    spreadData{r+row_offset, c+1} = data_str;
                end
            end


        end
    end
%%
    
    %% miscStats -- display for drifting & flashed
    if isfield(Sd, 'miscStats')
        dStats = Sd.miscStats;        
    else
        dStats = struct;    
    end
    dStatFields = fieldnames(dStats);                
    
    if isfield(Sf, 'miscStats')
        fStats = Sf.miscStats;         
    else
        fStats = struct;        
    end
    fStatFields = fieldnames(fStats);               
    
    statFields = uniqueInOrder([dStatFields; fStatFields]);
    
    
    spreadData{end+2,1} = 'MISC STATS FROM FILE';
    M = size(spreadData, 1);
    idx = 1;
    for si = 1:length(statFields)
        %%
        stat_i = statFields{si};
        spreadData{M+idx,1} = stat_i;
        haveD = isfield(dStats, stat_i);
        haveF = isfield(fStats, stat_i);
        
        
        if haveD            
            subFields = fieldnames( dStats.(stat_i) );
            nMaxD = max(structfun(@length, dStats.(stat_i)));
            if any(cellfun(@(fn) ~isempty(strfind(fn, '_p')), subFields))
                nMaxD = 2;
            end
            D_idx = 2; %:2+nMaxD-1;
        else
            nMaxD = 0;            
            subFields = fieldnames( fStats.(stat_i) );            
        end
        
        if haveF
            nMaxF = max(structfun(@length, fStats.(stat_i)));
            F_idx = 2+nMaxD; % + [0 : nMaxF-1];
        else
            nMaxF = 0;
        end
        
        
        D_header = iff(haveD, ['D', repmat({' '}, 1, nMaxD-1) ], {}); %D_idx = 2; %iff(haveD, 3, 0);
        F_header = iff(haveF, ['F', {}],                         {});% F_idx = 2+nMaxD;
        headerFields = [stat_i, D_header, F_header];
        [spreadData{M+idx, 1:length(headerFields)}] = headerFields{:};
        idx = idx+1;
        for sfi = 1:length(subFields)
            sub_fld_i = subFields{sfi};
            spreadData{M+idx, 1} = sub_fld_i;
            haveD_now = haveD && isfield(dStats.(stat_i), sub_fld_i);
            haveF_now = haveF && isfield(fStats.(stat_i), sub_fld_i);
            if haveD_now
                d_vals_str = toStr( dStats.(stat_i).(sub_fld_i), stat_i, sub_fld_i );
                for di = 1:length(d_vals_str)                                        
                    spreadData{M+idx, D_idx+di-1} = d_vals_str{di};
                end
            end
            if haveF_now
                f_vals_str = toStr( fStats.(stat_i).(sub_fld_i), stat_i, sub_fld_i );
                for fi = 1:length(f_vals_str)
                    spreadData{M+idx, F_idx+fi-1} = f_vals_str{fi};
                end                
            end                            
            
            idx = idx +1;
        end
        idx = idx +1;
        
%         spreadData{M+si,1} = statFields{i};
        
    end
        
    if opt.includeComments
        if isfield(Sd, 'comments')
            spreadData{end+1,1} = '';
            spreadData{end+1,1} = ' *-*-*-* Comments from drifting grating file *-*-*-*';
            for i = 1:length(Sd.comments)
                spreadData{end+1,2} = Sd.comments{i}; %#ok<AGROW>
            end
        end
        if isfield(Sf, 'comments')
            spreadData{end+1,1} = '';
            spreadData{end+1,1} = ' *-*-*-*Comments from flashed grating file *-*-*-*';
            for i = 1:length(Sf.comments)
                spreadData{end+1,2} = Sf.comments{i}; %#ok<AGROW>
            end
        end    
    end
    
    spreadData{end+1,1} = sprintf('Last updated: %s', datestr(now, 'mm/dd/yy HH:MM AM'));
    nBlankRowsToAdd = 20;
    [spreadData{end+[1:nBlankRowsToAdd], :}] = deal([]);
    
    timeWindow_str = iff(strcmp(timeWindow, 'best'), '_best', sprintf('_%d_%d', timeWindow));
    
    switch cmpType
        case 'phase',  
            sheetName = ['Ph_' statType timeWindow_str];
        case 'degree', 
            oe_mode_str = switchh(oe_mode, {'aa', 'oe_same', 'oe_diff'}, {'', '_OEs', '_OEd'});
            bccType_str = switchh(bccType, {'samePen', 'sameAnimal', 'diffPen-sameAnimal', 'full'}, ...
                                           {'_WP',      '_WA',        '_BP-WA',              ''});
                        
            sheetName = ['D_' statType iff(subtractSpont, '_SS', '_SI') bccType_str  iff(preserveSC, '_pSC', '') iff(preserveAB, '_pAB', ''), oe_mode_str ];
    end    
    
        warning('off', 'MATLAB:xlswrite:AddSheet');
    
    
    try
        xlswrite(excelFile, spreadData, sheetName);
    catch Merr
        if strcmp(Merr.identifier, 'MATLAB:xlswrite:LockedFile');
            beep;
            [~, file_nm, ext] = fileparts(excelFile);
            fprintf('Please close excel file:  "%s%s"  and try again ... \n', file_nm, ext);
            keyboard;
            xlswrite(excelFile, spreadData, sheetName);
        else
            rethrow(Merr)
        end
    end
    [~, excelFile_name, file_ext] = fileparts(excelFile);
    fprintf('  => Updated sheet "%s" in excel file %s%s.\n\n', sheetName, excelFile_name, file_ext);
    3;

end
    

function strC = toStr( vals, stat_i, fld_i)
    if iscellstr(vals) % column names
        strC = vals;
        return;
    end

%     if ~isempty(strfind(fld_i, '_p'))  % p-value
    if length(fld_i) > 1 && strcmp(fld_i(end-1:end), '_p')  % p-value
        assert(length(vals) == 1);
        if vals > 0.5
            strC = { sprintf('%.1f;', vals) };
        elseif vals > 0.1
            strC = { sprintf('%.2f;', vals) };
        elseif vals > 0.01
            strC = { sprintf('%.2f;', vals) };
        elseif vals > 0.001
            strC = { sprintf('%.3f;', vals) };
        else
            s_full = sprintf('%.0e', vals);
            [s1, s2] = strtok(s_full, 'e');
%             if vals > 0.00001 && ~strcmp(s1, '1')  % don't want 1x10                
%                 s1 = [s1 '×10'];
%             else
                s1 = ['10'];
%             end            
            s2 = strrep(s2(2:end), '0', '');            
            strC = { s1, s2 };            
        end
        return;
    end
    
            
    if ~iscell(vals)
        vals = num2cell(vals);
    end
        
        
    if ~isempty(strfind(fld_i, '_ci'))
        fmts = {'%.3f', '%.3f'};
        addColon = 1;

    elseif strncmp(fld_i, 'frac', 4)
        fmts = {'%.4f'};
        addColon = 1;

        
    elseif ~isempty(strfind(fld_i, 'n_')) || ~isempty(strfind(fld_i, 'num_'))
        fmts = repmat({'%d'}, 1, length(vals));
        addColon = 0;
        
    elseif strcmp(stat_i, 'outlierSigStats')
        fmts = {'%.3f', '%.3f'};
        addColon = 1;
        
    elseif ~isempty(strfind(fld_i, '_cc')) %% strncmp(stat_i, 'corr_ori', 7) ||  
        fmts = {'%.3f'}; % corr coeff, p_val
        addColon = 1;
        
    elseif strncmp(fld_i, 'spf_m', 5) || strcmp(fld_i, 'range_octaves')  % spf_min/spf_max
        fmts = {'%.3f'}; % corr coeff, p_val
        addColon = 1;        
        
    elseif ~isempty(strfind(fld_i, '_Wrcc'))         
        fmts = {'%.2f'}; % corr coeff, p_val
        addColon = 1;        
        
    else
        fmts = repmat({'%d'}, 1, length(vals));
        addColon = 0;
    end
        
    if addColon
        strC = cellfun(@(v, fmt) sprintf([fmt ';'], v), vals, fmts, 'un', 0);
    else
        strC = cellfun(@(v, fmt) sprintf([fmt ], v), vals, fmts, 'un', 0);
    end

end



function s = add_colon(s)
    if ~isnan(str2double(s))
        s = [s ';'];
    end
end

function [filename, opt] = getStatFileName(cmpType, gratingType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode, statType)

    switch statType
        case {'cellStats', 'cellDistribs', 'pairDiffs', 'pairStats'},      preserveSC = 0;
        case {'simpleComplexCells', 'simpleComplexPairs'}, preserveSC = 1;  %'scCellStats', 'scPairStats'
    end    
    oe_mode_type = switchh(cmpType, {'phase', 'degree'}, {'phaseOEmode', 'degreeOEmode'});
    opt = struct('cmpType', cmpType, 'subtractSpont', subtractSpont, 'gratingType', gratingType, 'preserveSimpleComplex', preserveSC, 'preserveAligned', preserveAB, 'bccType', bccType, 'timeWindow', timeWindow, oe_mode_type, oe_mode);
    filename = getFileName(statType, [], [], opt); 
    
end

function confirmCorrectFile(Sf, opt)    
    assert(strcmp(Sf.cmpType, opt.cmpType));
    assert(strcmp(Sf.gratingType, opt.gratingType));
    if isfield(Sf, 'subtractSpont')
        assert(isequal(Sf.subtractSpont, opt.subtractSpont));
    end
    if isfield(Sf, 'preserveSimpleComplex')
        assert(isequal(Sf.preserveSimpleComplex, opt.preserveSimpleComplex));
    end
    if isfield(Sf, 'bccType')
        assert(isequal(Sf.bccType, opt.bccType));
    end
end

function tf = excludeField(measureName, gratingType)
    tf = 0;
    switch gratingType
        case 'flashed',            
            tf  = 0 && any( cellfun(@(s) ~isempty(strfind(measureName, s)), {'R_spont', 'R90_stim', 'R90_total'}));
            
        case 'drifting',

    end
end


function data_str = format_data(data_val, measureName, column_name, statType, opt)
    
    SC_pval_th = 0.001;
    isSCpval = strncmpi(column_name, 'pval_', 5) || ((length(column_name) >= 2) && strcmp(column_name(end-1:end), '_p'));
    isNPair = ~isempty(strfind(column_name, 'N1'));
    isCCMeasure = ~isempty(strfind(measureName, 'cc'));
    isOriOrDirMeasure = any( cellfun(@(s) ~isempty(strfind(measureName, s)), {'_ori', '_dir'}));
    isDphiMeasure = ~isempty(strfind(lower(measureName), 'dph'));
    isMCProbMeasure = any( cellfun(@(s) ~isempty(strfind(lower(column_name), s)), {'prob'}));
    isKSstat = any( cellfun(@(s) ~isempty(strfind(lower(column_name), s)), {'ksstat'}));
    isStandardPval = any( cellfun(@(s) ~isempty(strfind(lower(column_name), s)), {'ks', 'pval'})) || (length(column_name) > 1 && strcmp(column_name(end-1:end), '_p'));
    isRatio = ~isempty(strfind(lower(column_name), 'ratio'));
    isFrac = strncmp(column_name, 'frac', 4);
    
    
    doConvertSmallProb = 0; 
    doRemoveExtraZeros = 0;
    doRemovePrFromBcc = 1;
    addDegreeSign = 0;
    Ntrunc = [];
    
    if isSCpval
        if (data_val > 0) && (data_val < SC_pval_th)
            data_val = round(log10(data_val));
            fmt_code = '%d';
        else            
            fmt_code = '%.4f'; Ntrunc = 1;
        end
    elseif ischar(data_val)
        fmt_code = '%s';
    elseif isFrac
        fmt_code = '%.4f';
    elseif strcmp(column_name, 'N') || strncmpi(column_name, 'N_', 2)
        fmt_code = '%d';
    elseif isKSstat
        fmt_code = '%.2f';
    elseif isMCProbMeasure
        fmt_code = '%.4f'; Ntrunc = 1;
    elseif isCCMeasure
        fmt_code = '%.2f'; Ntrunc = []; 
    elseif (isOriOrDirMeasure && ~strcmp(statType, 'pairStats') && ~isRatio) ...  % ori & dir: only add in cellStats, pairDiffs, not pairStats (don't want degree sign with medianRatios, medianProbs, etc.)
            || (isDphiMeasure && ~isStandardPval)
        fmt_code = '%.2f'; Ntrunc = 2;
        addDegreeSign = 1;
        
    elseif isStandardPval
        fmt_code = '%.1g'; doRemoveExtraZeros = 1;
    else
        fmt_code = '%f'; Ntrunc = 2;
    end

    data_str = sprintf(fmt_code, data_val);
    
    tf_converted = 0;
    if isMCProbMeasure && doConvertSmallProb && opt.separateCellsForFlashedAndDrifting
        [data_str, tf_converted] = convertSmallProb(data_str);
    end
    if doRemoveExtraZeros
        data_str = removeExtraZerosFromExpStr(data_str);
    end
    if isNPair && doRemovePrFromBcc
        s = sscanf(data_str, '%d Pr');
        if s > 10000
            data_str = strrep(data_str, ' Pr', '');
        end
        3;
        
    end
    if ~isempty(Ntrunc) && ~tf_converted
        data_str = truncateDecimalAfterNSigDigits(data_str, Ntrunc);
    end
    if isMCProbMeasure && ~doConvertSmallProb 
        if strcmp(data_str, '0.0')
            data_str = '0';
        end
        
    end
    if addDegreeSign
        data_str = [data_str '°'];
    end
    
%     if strncmp(data_str, '0.0000', 6)
%         3;
%     end



end


function [str, tf_converted] = convertSmallProb(str)
    val = str2double(str);
    tf_converted = val < 1e-4;
    if tf_converted;
        str = '<1e-4';
    end
end


function fn_str = extractFilename(fn_withPath)
    [~, fn, ext] = fileparts(fn_withPath);
    fn_str = [fn ext];
end

