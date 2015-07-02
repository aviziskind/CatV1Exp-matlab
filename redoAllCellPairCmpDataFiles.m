function redoAllCellPairCmpDataFiles

    gratingType_vals = {'flashed', 'drifting'};
    
    redoAll = 0;
    redoOldFiles = 1;
        redoDate_cellFiles    = 735995.800658;    %sprintf('%.6f', now)
        redoDate_pairCmpFiles = 735995.800658;
    
    redoPairCmpFilesIfOlderThanCellFiles = 1;
       
    haveFile = @(filename, redoBeforeDate) exist(filename, 'file') && (~redoOldFiles || ~fileOlderThan(filename, redoBeforeDate)) && ~redoAll;
        
    haveFiles = @(fileList, redoBeforeDate) cellfun(@(fn) haveFile(fn, redoBeforeDate), fileList);

    doDegree = 1;
    doPhase = 0;
    
            
                         % Cmp          subtSpont  bccType   preserveSC  preserveAB timeWindow  odd/even   
    allCombinations = { ...
                         'degree',        1        'full'    0          0         'best',   'aa'   ;  ... standard degree of tuning test, SS.
                      ...'degree',        1        'full'    0          0         'stimw',   'aa'   ;  ... degree of tuning test, SS, comparing odd/even.
                         'degree',        0        'full'    0          0         'best',   'aa'   ;  ... standard degree of tuning test, SI.
                         'degree',        1        'samePen',0          0         'best',   'aa'   ;  ... only same penetration, SS (test for spf)
                         'degree',        0        'samePen',0          0         'best',   'aa'   ;  ... only same penetration, SI (test for spf)
                         'degree',        1        'full'    1          0         'best',   'aa'   ;  ... for simple/complex clustering analysis (SS)
                         'degree',        1        'full'    0          1         'best',   'aa'   ;  ... for test of clustering of preferred direction (SS)
                      ...'degree',        1        'samePen' 1          0         'best',   'aa'     ;  ... for simple/complex clustering of spatial frequency (SS)
                         'degree',        1        'BP-WA',  0          0         'best',   'aa'     ;  ... for simple/complex clustering of spatial frequency (SS)
                         'degree',        1        'full'    0          0         'best',   'oe_diff'   ;  ... degree of tuning test, SS, comparing odd/even (diff)
                         'degree',        1        'full'    0          0         'best',   'oe_same'   ;  ... degree of tuning test, SS, comparing odd/even (same)
                         'degree',        1        'full'    1          0         'best',   'oe_diff'   ;  ... degree of tuning test, SS, comparing odd/even (diff) - for SS pairs
                         'degree',        1        'full'    1          0         'best',   'oe_same'   ;  ... degree of tuning test, SS, comparing odd/even (same) - for SS pairs
                         

                         'phase',         0        'full'    0          0         'best',   'oe_diff'     ;  ... standard phase tuning. (no between site, and no Wrcc, so preserveSC and bccType are irrelevant)                   
                         'phase',         0        'full'    0          0         [29 62],  'oe_diff'    ;  ... phase tuning with fixed window
                         'phase',         0        'full'    0          0         [58 91],  'oe_diff'    ;  ... phase tuning with fixed window
                         'phase',         0        'full'    0          0         'best',   'oe_same'     ;  ... standard phase tuning, odd-odd, even-even
                         'phase',         0        'full'    0          0         'best',   'aa'          ;  ... standard phase tuning, all-all
                         }; 
                     
                     
     nCombinations = size(allCombinations, 1);
     for comb_i = 1:nCombinations
         [cmpType, subtractSpont, bccType, preserveSC, preserveAB, timeWindow, oe_mode] = allCombinations{comb_i, :};
         if strcmp(cmpType, 'degree') && ~doDegree
             continue;
         end
         if strcmp(cmpType, 'phase') && ~doPhase
             continue;
         end
         
         curCmpType(cmpType);
         curSubtractSpont(subtractSpont);
         curBccType(bccType);
         curPreserveSimpleComplex(preserveSC);
         curPreserveAligned(preserveAB);         
         curTimeWindow(timeWindow);
         
         if strcmp(cmpType, 'degree')
             curDegreeOEmode(oe_mode);
         end
         if strcmp(cmpType, 'phase')
             curPhaseOEmode(oe_mode);
         end

         window_str = iff(strcmp(timeWindow, 'best'), 'best window', sprintf('window:[%d-%d]', timeWindow));                 
         
         switch cmpType
             case 'degree', pairTypes = {'Wcc', 'Wrcc', 'Bcc'};
             case 'phase',  pairTypes = {'Wcc', 'Wscc'};
         end
         curPairTypes(pairTypes);
         
         for gt_i = 1:length(gratingType_vals)
             curGratingType(gratingType_vals{gt_i});
             if strcmp(gratingType_vals{gt_i}, 'drifting') && ~strcmp(timeWindow, 'best')
                 continue; % window option only applies for flashed grating. 
             end
             
                 
             fprintf('*************************************************************************************************************\n');
             fprintf('* Set %d/%d. CmpType: %s. SubtSpont: %d. bccType: %s. preserveSC: %d. preserveAB: %d. Grating: %s. OE: %s. [window:%s]\n', ...
                 comb_i, nCombinations, cmpType, subtractSpont, bccType, preserveSC, preserveAB, gratingType_vals{gt_i}, oe_mode, window_str);
             fprintf('*************************************************************************************************************\n');

             %%
             
             cells_files = getFileNamesNeeded('osps', cmpType, preserveAB);
             if any(~haveFiles(cells_files, redoDate_cellFiles));
                generateGratingCellsDatafile;
             end                     

             pair_files = getFileNamesNeeded('pairs', cmpType, preserveAB);
             if any(~haveFiles(pair_files, redoDate_pairCmpFiles)) || (any(cellfun(@(fn) fileOlderThan(fn, cells_files{1}), pair_files)) && redoPairCmpFilesIfOlderThanCellFiles)
                 generateGratingPairsDatafile;
             end
             
             cmp_files = getFileNamesNeeded('comparisons', cmpType, preserveAB);
             if any(~haveFiles(cmp_files, redoDate_pairCmpFiles)) || (any(cellfun(@(fn) fileOlderThan(fn, pair_files{1}), cmp_files)) && redoPairCmpFilesIfOlderThanCellFiles)
                 generateGratingComparisonsDatafile;
             end             

         end
     end


end

function allFileNames = getFileNamesNeeded(fileType, cmpType, preserveAB)

    switch cmpType
        case 'degree',            
              switch fileType                          
                  case 'osps',      
                     ori_cells_file = getFileName('osps', '_ori');
                     spf_cells_file = getFileName('osps', '_spf');
                     all_cells_file = getFileName('osps', '_all');
                     allFileNames = {ori_cells_file, spf_cells_file, all_cells_file};
                      
            
                case 'pairs',            
                    ori_pairs_file = getFileName('pairs', '_ori');
                    spf_pairs_file = getFileName('pairs', '_spf');
                    allFileNames = {ori_pairs_file, spf_pairs_file};
                    if preserveAB
                        allFileNames = allFileNames(1);
                    end
                    
                case 'comparisons',
                    ori_cmp_file = getFileName('comparisons', '_ori');
                    spf_cmp_file = getFileName('comparisons', '_spf');
                    allFileNames = {ori_cmp_file, spf_cmp_file};
                    if preserveAB
                        allFileNames = allFileNames(1);
                    end

              end
                                
            
        case 'phase',
            
            switch fileType                          
                  case 'osps',      
                     cells_file = getFileName('osps');
                     allFileNames = {cells_file};                      
            
                case 'pairs',            
                     pairs_file = getFileName('pairs');
                     allFileNames = {pairs_file};
                    
                case 'comparisons',
                    cmp_file = getFileName('comparisons');
                    allFileNames = {cmp_file};
            end
    end

end

             
           


%     drifting & flashed:
%         degree of tuning:        
%             subtractSpont = 0;
%                 preserveSimpleComplex = 0;
%                     bccType = 'full'
%                     bccType = 'Pen'
% 
%                 preserveSimpleComplex = 1; % just for S/C statistics pairs.
%                     bccType = 'full'
%                     
%             subtractSpont = 1;
%                 preserveSimpleComplex = 0;
%                     bccType = 'full'
%                     
%         phase tuning:        
%             preserveSimpleComplex = 0;
%                 subtractSpont = 0;
%                     bccType = 'full'

    