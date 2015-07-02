function [sorting_str_out, cmp_str_out] = curStatus

    % whether are using the spikesorting saved in the database (manual sorting done using Spiker)
    matchDB = curMatchDB;
    
    % which clustering level - 'clusters', 'clustersPruned', 'or 'cells'
    clustGroupType = curGroupingType('');
    
    sortingFeatures = curSortingFeatures('');
    
    if matchDB
        sortingFeatures = '(Manual)';
    end
        
    gratingType = curGratingType(''); % flashed or drifting
    
    cmpType = curCmpType(''); % phase tuning or degree of tuning
    cmpType_str = switchh(cmpType, {'phase', 'degree', 'psth'}, {'phase tuning', 'degree of tuning', 'psth similarity'});
        
%     simpComp = curSimpleComplex('');
%     simpComp = 'all';
%     simpComp_str = switchh(simpComp, {'simple', 'complex', 'all'}, {'Simple Cells only', 'Complex Cells only', '(Simple and Complex cells)'});
    
%     cmpType = curCmpType; % phase tuning comparisons (usually) or psth comparisons
    
    pt = curPairTypes(''); % within site only, between site only, or both.
    pairTypes = cellstr2csslist(pt, ', ');
    
    sorting_str = sprintf('[Match DB: %s]  Using "%s" clustering. Sorting features: "%s"', iff(matchDB, 'YES', 'NO'), clustGroupType, sortingFeatures);
%     cmp_str = sprintf('Comparing %s for *%s* gratings (%s), using %s pairs', cmpType_str, gratingType, simpComp_str, pairTypes);
    cmp_str = sprintf('Comparing %s for *%s* gratings, using %s pairs', cmpType_str, gratingType, pairTypes);
    
    if nargout == 0
        fprintf('%s\n%s\n', sorting_str, cmp_str);
    else
        sorting_str_out = sorting_str;
        cmp_str_out = cmp_str;
    end
    
    % note : these 'cur-' files are 'short cut' files
    % curCellsType
    % curCellGroupFile

end