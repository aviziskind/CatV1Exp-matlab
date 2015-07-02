function exploreComparisonOfGratingCells

    printProgress = false;
    useMeasurePvals = false;
    legendFontSize = 13;
    
    showControlPanel       = 1;
    makeParameterPlots     = 0;
    makeBarGraphOfResults  = 0;

    skipBigHeading = true;
    skipYlabels = true;
    showLocLabelsIfOneLabel = false;
    showSigTestsOnPlots = true;

    hideLegend = true;
        fastLegendUpdates = false;    
    
    [gratingType, gratingType_s] = curGratingType;  % FLASHED_GRATING = 1; DRIFTING_GRATING = 2;
    curCmpType('phase');
    [cmpType, cmpType_s] = curCmpType;    
    [pt_ids, pairTypes] = curPairTypes;    
        
    ospDatafile = getFileName('osps');
    pairDatafile= getFileName('pairs');
    cmpDatafile = getFileName('comparisons');
    
    fprintf('Loading ... '); tic;    
    S1 = load(ospDatafile);        
    allCells = S1.allCells;
    clear S1;
    nUnits = length(allCells);    
    
    S2 = load(pairDatafile);    
    [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs, Wrcc_pairIdxs, Wrcm_pairIdxs] = ...
        deal(S2.Wcc_idxs, S2.Wcm_idxs, S2.Bcc_idxs, S2.Bcm_idxs, S2.Wrcc_idxs, S2.Wrcm_idxs);   
    clear S2;
    
    S3 = load(cmpDatafile); 
    [pairData, S, pairTypes, measures, locations, opt] = deal(S3.pairData, S3.allStatsC, S3.pairTypes, S3.measures, S3.locations, S3.opt);
    clear S3;
%     flds_pd = fieldnames(pairData);
%     for fld_i = 1:length(flds_pd)
%         pairData_cat.(flds_pd{fld_i}) = vertcat( pairData.( flds_pd{fld_i} ) );
%     end
    
    
    fprintf('done. '); toc;
        
%     pairTypes = {'Wcc'};
%     pairTypes = {'Wcc', 'Wrcc', 'Wcm', 'Wrcm'};    
    
    removeWrcc = 1;
    removeWscc = 1;
    if removeWrcc
        pairTypes(find(strcmp(pairTypes, 'Wrcc'), 1)) = [];
    end
    if removeWscc
        pairTypes(find(strcmp(pairTypes, 'Wscc'), 1)) = [];
    end
    

    % Pair Types
    
    allPairIdxs = {Wcc_pairIdxs, Wrcc_pairIdxs, Bcc_pairIdxs, Wcm_pairIdxs, Wrcm_pairIdxs, Bcm_pairIdxs};
    allPairIdxs_list = {Wcc_pairIdxs, uniqueInts(Wrcc_pairIdxs), Bcc_pairIdxs, Wcm_pairIdxs, uniqueInts(Wrcm_pairIdxs), Bcm_pairIdxs};
    allPairTypes = {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
    pairTypes = pairTypes(ord(cellfun(@(s) find(strcmp(s, allPairTypes)), pairTypes)));
    
%     allPairTypeLabels = {'Within site: cell-cell', 'Within-site: cell-cell (rand)',  'Between site: cell-cell', ...
%                          'Within site: cell-multiunit', 'Within site: cell-multiunit (rand)', 'Between site: cell-multiunit'};
    allPairTypeLabels = {'', 'Within-site: cell-cell (rand)',  '(Control)', ...
                         'Within site: cell-multiunit', 'Within site: cell-multiunit (rand)', 'Between site: cell-multiunit'};
    pairTypesAvailable = cellfun(@(pr) any ( strcmp(pairTypes, pr)), allPairTypes);
    pairTypeLabels = cellfun( @(pr) allPairTypeLabels{find(strcmp(allPairTypes, pr))},   pairTypes, 'un', false);
    pairIdxs = cellfun( @(pr) allPairIdxs{find(strcmp(allPairTypes, pr))},  pairTypes, 'un', false);

    rescaleDphiDistribution = true;    
    rescaleBccDistribution = true;
    
    rescaleBccDistribution = rescaleBccDistribution && any( strcmp(pairTypes, 'Wcc') ) && any( strcmp(pairTypes, 'Bcc') );
    
    nP = length(pairTypes);
    pairTypesActive = true(1,nP);
    pairIdxsResampled = find(strcmp(cellfun(@class, pairIdxs, 'un', false), 'cell'));
    if ~isempty(pairIdxsResampled)
        nMaxResamples = length(pairIdxs{pairIdxsResampled(1)});
    else
        nMaxResamples = 1;
    end
    nResamples = min(10, nMaxResamples);
    red_col = [1, 0, 0];  pink_col = [1, 0, .5];   purple_col = [1 0 1];
    green_col = [0 1 0];  aqua_col = [0 1 .5];     turq_col = [0 1 1];
    allPairColors = {red_col, pink_col, purple_col,   green_col, aqua_col, turq_col};
    pairColors = allPairColors(pairTypesAvailable);
    
    pairIdxList = uniqueInts(allPairIdxs_list(pairTypesAvailable));
    nPairs = length(pairIdxList);
    
    % Locations     
    allLocations      = {'maxR1',  'maxR2',  'mean12',  'maxR1xR2', 'maxMinFracR', 'maxMU', 'wgtedSum_all', 'wgtedSum_top10',  'wgtedSum_top20', 'wgtedSum_above90', 'wgtedSum_above50'};
    allLocationLabels = {'max R1', 'max R2', 'mean 1&2', 'max R1xR2', 'max min(R1,R2)', 'max R_{MU}', '[all]', '[top10]', '[top20]', '[above90%]', '[above50%]'}; % change back to "all"    
    
    extraLocations = setdiff(locations, allLocations);
    allLocations = [allLocations, extraLocations];
    allLocationLabels = [allLocationLabels, extraLocations];
    orig_loc_idx = cellfun(@(s) find(strcmp(s, allLocations)), locations);
    assert(isequal(locations, allLocations(orig_loc_idx))); 
    locationLabels = allLocationLabels(orig_loc_idx);
%     locationsAvailable = cellfun(@(loc) any ( strcmp(locations, loc)), allLocations);                        
    
    nLoc = length(locations);    
    locationsActive = true(1,nLoc);

    % Measures
    allMeasures      = {'dot', 'cc',  'rho',  'tau',  'dphi',      'dF1',     'STA_cc',    'MID_cc',   'MID_ovlp_cc',   'MID_fit_cc', 'dPh_rel' };
    allMeasureLabels = {'dot', 'cc', '\rho', '\tau', '\Delta\phi', '\DeltaF1', 'cc_{STA}', 'cc_{MID}', 'cc_{MID-ovlp}', 'cc_{MID-fit}', 'dPh_{rel}'};
    allMeasureNullMeans = [.5,  0,   0,       0,         90,          90,      0,           0,         0,                0,            90];
    measure_neword = ord(cellfun(@(s) find(strcmp(s, allMeasures)), measures));
    measures = measures(measure_neword);
    S = S(:,measure_neword);
    
%     measuresAvailable = cellfun(@(ms) any ( strcmp(measures, ms)), allMeasures);
    measureInds = cellfun( @(ms) find(strcmp(allMeasures, ms)),   measures);           
    measureLabels = allMeasureLabels(measureInds);
    mMullMeans = allMeasureNullMeans(measureInds);
    nM = length(measures);
    measuresActive = true(1,nM);    
    idx_dphi = find(strcmp('dphi', measures));
    idx_dF1 = find(strcmp('dF1', measures));
    if ~isempty(idx_dphi)
        dPhi_binE = S{1, idx_dphi}.binEdges;
    end
    if ~isempty(idx_dF1)
        dF1_binE = S{1, idx_dF1}.binEdges;
    end
%     nF1_bins = length(S{1,find(strcmp(measures, 'dF1'))}.binEdges)-1;

%     nPh_sensitiveMeasures = {'cc', 'rho', 'tau', 'dphi'};
    
    % Filters / Thresholds / Limits    
    fontsize = 13;
    nMaxSigTests = 7;
    fmt = '%.2g';
%     str_mn_std  = @(x,w) [num2str(mean_wgt(x, w), fmt) ' \pm ' num2str(stderr_wgt(x, w), fmt) ' (' num2str( nnz(~isnan(x))) ')' ];
    Nstr = @(n) iff(isempty(n), '', iff(round(n) == n, sprintf(' (N = %d)', n), sprintf(' (%.1f)', n) ) );    % 
    str_mn_std_n = @(mn,st, n) sprintf( [fmt ' \\pm ' fmt ' %s'], mn, st, iff(isempty(n), '', Nstr(n))); % skip median
%     str_mn_std2 = @(mn,st, md, n) sprintf( [fmt ' \\pm ' fmt '; ' fmt '%s'], mn, st, md, dispN(n));  % include median 
%     str_mn_std2 = @(mn,st, n) [sprintf( [fmt ' \\pm ' fmt], mn, st), iff(~isempty(n), '', '' )];
    str_pval    = @(p)  sprintf('p = %.2g', p);
    makePvalStr = @(pvals, pval_labels) cellstr2csslist( cellfun(@(p,s) [s ' = ' num2str(p, '%.2g')], num2cell(pvals(1:length(pval_labels))), pval_labels, 'un', 0) );
    wrp = @(x) [x(1:end); x(1)];
    ext = @(x) [x(1:end); x(end)+diff(x(1:2))];
    shft = @(x,n) x([end-n+1:end, 1:end-n]);

%{    
    function y = locCmp(val, cmpMode)
        if strncmp(cmpMode, 'same', 4)
            y = val(:,1);
        elseif strncmp(cmpMode, 'diff', 4)
            y = ~val(:,1);
        elseif strncmp(cmpMode, 'auto', 4)
            y = val(:,2);
        end        
    end
%}                                


%     title( str_mn_std (cc(:,x_i)), 'fontsize', fontsize-1 );
%     th      = struct('minF1oDC', .5,  'minFracOfMaxes', .5, 'pval', 1,  'min_oriSelPval', .05, 'min_oriRepPval', .05, 'min_spfRepPval', .05);
%     th_prev = struct('minF1oDC', nan,  'minFracOfMaxes', nan, 'pval', nan,  'min_oriSelPval', nan,  'min_oriRepPval', nan,  'min_spfRepPval', nan);              
%     alpha_05 = -log10(.05);
    allNPhases = double( unique(pairData.n_phases) );    
%     animalCmpStrs = {'same', 'diff', 'auto'};    
    
    anyPvals = cellfun(@(s) isfield(s, 'pval'), S);
    if ~any(anyPvals(:)), useMeasurePvals = false; end

    rng01 = 0:.01:1;
    rng02 = 0:.01:2;
    rng_1_1 = -1:.01:1;
    locCmpStrs = {'same', 'diff'};
%     SCtypesStr = {'S-S', 'S-C', 'C-C'};
    SCtypesStr = {'S-S', 'S-C', 'C-C'};
    pvals = [0:.01:1];
    pvalsRng = iff(useMeasurePvals, pvals, 0);
    logPvalRng = [0:.1:10];
    allNPhases_S = arrayfun(@num2str, allNPhases, 'un', 0);
    
%     pairFilterNames     = {'min_rep_cc_p', 'n_phases',    'SCtype_pref',  'SCtype_cmp', 'minF1oDC_pref_lo', 'minF1oDC_pref_hi', 'maxF1oDC_pref_lo', 'maxF1oDC_pref_hi', 'minJackCC_lo', 'minJackCC_hi', 'minRsqr_lo', 'minRsqr_hi',  'animalCmp', 'penetrCmp', 'locCmp'   };  % full name - goes on var label boxes
%     pairFilterVarNames  = {'min_rep_cc_p', 'n_phases',    'SCtype_pref',  'SCtype_cmp', 'minF1oDC_pref',    'minF1oDC_pref',    'maxF1oDC_pref',    'maxF1oDC_pref',    'min_jackCC',   'min_jackCC',   'min_rqr',    'min_rsqr',    'animalCmp', 'penetrCmp', 'locCmp'   };  % name of actual variable that is affected
%     pairFilterVarOp     = {@ge,            @eq,           @eq,             @eq,          @ge,                @le,                @ge,                @le,               @ge,             @le,           @ge,           @le,          @eq,         @eq,          @eq       };     
%     pairFilterVarInit   = { 0,              allNPhases(1), 0,              0,           .75,                 1,                 .5,                 1,                  0.4,             1,              0.4,          1,             1,           1,            1         };
%     pairFilterRange     = { logPvalRng,     allNPhases_S,  SCtypesStr,     SCtypesStr,   rng02,              rng02,              rng02,             rng02,              rng_1_1,         rng_1_1,        rng01,        rng01,         locCmpStrs,  locCmpStrs,   locCmpStrs};
%     pairFilterGroup     = { 1,              1,             2,              2,            3,                  3,                  4,                 4,                  5,               5,              5,            5,             6,           6,            6         };
%     pairFilterLabel     = { 'rep_p_av',    'n_phases',    'S/C_pref',     'S/C_cmp',    'mnF1PrefLo',       'mnF1PrefHi',       'mxF1PrefLo',      'mxF1PrefHi',       'minJackCC_lo',  'minJackCC_hi',  'minRsqr_lo', 'minRsqr_hi', 'Animal',    'Penetration','Location' };  % short name - goes on vector labels.


%     locFilterNames    = {'minF1oDC_cmp_lo', 'minF1oDC_cmp_hi', 'maxF1oDC_cmp_lo', 'maxF1oDC_cmp_hi', 'minFracR_lo',    'minFracR_hi',    'minSumPhs',  'maxSumPhs', 'minPhTC_cc',    'minPhTC_dot'   , 'cc_atMaxDphi' };
%     locFilterVarNames = {'minF1oDC_cmp',    'minF1oDC_cmp',    'maxF1oDC_cmp',    'maxF1oDC_cmp',    'minFracOfMaxes', 'minFracOfMaxes', 'minSumPhs',  'minSumPhs', 'minPhaseTC_cc', 'minPhaseTC_dot', 'cc_atMaxDphi'};    
%     locFilterVarOp    = {@ge ,             @le,                @ge ,              @le,               @ge,             @le,               @ge,         @le,             @ge,             @ge             , @ge        };
%     locFilterVarInit  = {.75,               1,                 .5,                .7,                .65,              1,                 1,           100            -1,              -1                 0          };
%     locFilterVarRange = { rng02,            rng02,              rng02,             rng02,             rng01,           rng01,            [0:1:360],   [0:1:360]      rng_1_1,          rng_1_1           rng_1_1     };
%     locFilterVarGroup = { 1,                1,                  2,                 2,                 3,               3,                 4,            4,             4,                4              ,  4         };
%     locFilterLabel    = {'mnF1CmpLo',      'mnF1CmpHi',        'mxF1CmpLo',       'mxF1CmpHi',       'minFracR_lo',   'minFracR_hi',     'mnSmPhs',     'mxSmPhs'   'phTC_cc',        'phTC_dot'      ,  'dp_cc',    }; 
        
    
    filters_fields = ...
        {'filterName',        'varName',        'op',   'init',         'range',     'group', 'label',    'type'};
    
    % filters on pair-dependent variables
    pairFilters_Cell = {...
    ... Filter Name           VarName            Op      InitValue      VarRange        Group Label
        'min_rep_cc_p'        'min_rep_cc_p'     @ge,    0,             logPvalRng      1,    'rep_p_av'    
        'n_phases'            'n_phases'         @eq,    allNPhases(1), allNPhases_S    1,    'n_phases'    
        'SCtype_pref'         'SCtype_pref'      @eq,    0              SCtypesStr      2,    'S/C_pref'    
        'SCtype_cmp'          'SCtype_cmp'       @eq,    0              SCtypesStr      2,    'S/C_cmp'     
        'minF1oDC_pref_lo'    'minF1oDC_pref'    @ge,    0.75           rng02           3,    'mnF1PrefLo'  
        'minF1oDC_pref_hi'    'minF1oDC_pref'    @le,    1              rng02           3,    'mnF1PrefHi'  
        'maxF1oDC_pref_lo'    'maxF1oDC_pref'    @ge,    0.5            rng02           4,    'mxF1PrefLo'  
        'maxF1oDC_pref_hi'    'maxF1oDC_pref'    @le,    1              rng02           4,    'mxF1PrefHi'  
        'JackCC_lo'           'min_jackCC'       @ge,    0.4            rng_1_1         5,    'minJackCC_lo'
        'JackCC_hi'           'max_jackCC'       @le,    1              rng_1_1         5,    'minJackCC_hi'
        'Rsqr_lo'             'min_rsqr'         @ge,    0.4            rng01           5,   'minRsqr_lo'  
        'Rsqr_hi'             'max_rsqr'         @le,    1              rng01           5,    'minRsqr_hi'  
        'animalCmp'           'animalCmp'        @eq,    1              locCmpStrs      6,    'Animal'      
        'penetrCmp'           'penetrCmp'        @eq,    1              locCmpStrs      6,    'Penetration' 
        'locCmp'              'locCmp'           @eq,    1              locCmpStrs      6,    'Location'
        'minOSI'              'min_OSI'          @ge,    .1             rng01           7,    'min_OSI'      
        'maxOSI'              'max_OSI'          @le,    .9             rng01           7,    'max_OSI'      
        'negAmpOvlp',         'negAmps_overlap', @ge,    2              [0:1:50],       7,    'Overlap',...
        };

        pairFilters_Cell(:, end+1) = {'pair'};
        
    pairFilters = cell2struct(pairFilters_Cell, filters_fields, 2);
    
    % filters on location-dependent variables
    locFilters_Cell = {...    
    ... Filter Name           VarName            Op   InitValue    VarRange   Group   Label
        'minF1oDC_cmp_lo'     'minF1oDC_cmp'     @ge     0.75      rng02        1    'mnF1CmpLo'  
        'minF1oDC_cmp_hi'     'minF1oDC_cmp'     @le     1         rng02        1    'mnF1CmpHi'  
        'maxF1oDC_cmp_lo'     'maxF1oDC_cmp'     @ge     0.5       rng02        2    'mxF1CmpLo'  
        'maxF1oDC_cmp_hi'     'maxF1oDC_cmp'     @le     0.7       rng02        2    'mxF1CmpHi'  
        'minFracR_lo'         'minFracOfMaxes'   @ge     0.65      rng01        3    'minFracR_lo'
        'minFracR_hi'         'minFracOfMaxes'   @le     1         rng01        3    'minFracR_hi'
        'minSumPhs'           'minSumPhs'        @ge     1         [0:1:360]    4    'mnSmPhs'    
        'maxSumPhs'           'minSumPhs'        @le     100       [0:1:360]    4    'mxSmPhs'    
        'minPhTC_cc'          'minPhaseTC_cc'    @ge     -1        rng_1_1      4    'phTC_cc'    
        'minPhTC_dot'         'minPhaseTC_dot'   @ge     -1        rng_1_1      4    'phTC_dot'   
        'cc_atMaxDphi'        'cc_atMaxDphi'     @ge     0         rng_1_1      4    'dp_cc'};    
        locFilters_Cell(:, end+1) = {'loc'};    
    locFilters = cell2struct(locFilters_Cell, filters_fields, 2);

    % filters on the actual value
    valFilters_Cell = {...
        'pval',             'pval',              @le,    1     pvalsRng,    1   'pval' };
        valFilters_Cell(:, end+1) = {'val'};    
    
    valFilters = cell2struct(valFilters_Cell, filters_fields, 2);

    allFilters = [pairFilters; locFilters; valFilters];
    
    filterTypes  = {allFilters.type};
    filterFields = {allFilters.filterName};% {pairFilterVars{:}, locFilterVars{:}, valFilterVars{:}};
    filterOps    = {allFilters.op}; % {pairFilterVarOp{:}, locFilterVarOp{:}, valFilterVarOp{:}};

    pairFilterInds = find(strcmp('pair', filterTypes));
    locFilterInds  = find(strcmp('loc',  filterTypes));
    valFilterInds  = find(strcmp('val',  filterTypes));    
    
    filterActive = false(1, length(allFilters));
    setFilterActive_prev = cell(nP, nLoc, nM);
    setFilterActive_prev(:) = {false(1, length(allFilters))};
        
    initThVals     = {allFilters.init}; %[pairFilterVarInit, locFilterVarInit, valFilterVarInit];
    th_fields      = [filterFields; initThVals];
    th_prev_fields = [filterFields; num2cell(nan(size(initThVals), class(initThVals{1})))];
    th = struct(th_fields{:});
    setTh_prev = cell(nP, nLoc, nM);
    setTh_prev(:) = {struct(th_prev_fields{:})}; 
        
    nResamplesPrev = nResamples;    
%     pairTypesActive_prev = [];
    
%     filterMask_P = {}; 
%     filterMasks_PL = cell(nLoc, 1);    
    
    pairSubMasks = struct;    
    locSubMasks = cell(1, nLoc);        
    valSubMasks = cell(nLoc, nM);
    
    pairMask = true(nPairs, 1);    
    locMasks(1:nLoc) = {true(nPairs, 1)};
    valMasks(1:nLoc, 1:nM) = {true(nPairs, 1)};
    
    dataMasks = cell(nLoc, nM);
    
    % Color coding 
    pvalLims = [0, .01, .05, 1];
    F1oDC_thresholds = [0, .5, .75, 1, 2.001];
    maxMinFrac_thresholds = [0, .25, .5, .75, 1];
    
    hasMU = ~cellfun(@isempty, strfind(pairTypes, 'm'));
%                 SCtypes = [2, 0, 1.5, 1]; % [0, 1, 1.5, 2] = [cC-mC / cC-mS / cS-mC / cS-mS]

%     F1oDCpairingLabels(~hasMU) = {{'S-S', 'C-C', '', 'S-C',}};  % for F1/DC pairing type (without multi-unit)    
    F1oDCpairingLabels(~hasMU) = {{'Simple-Simple_{pref}', 'Complex-Complex_{pref}', '', 'Simple-Complex_{pref}'}};  % for F1/DC pairing type (without multi-unit)    
    F1oDCpairingLabels(hasMU)  = {{'mS-cS', 'mC-cC', 'mC-cS', 'mS-cC'}};                       % for F1/DC pairing type (with multi-unit)
    F1oDCthresholdLabels(1:nP) = {legendarray('min F1/DC > ', F1oDC_thresholds(1:end-1) )};  % for F1oDC_threshold
    MinFracThresholdLabels(1:nP) = {legendarray('min Frac > ', maxMinFrac_thresholds(1:end-1) )};           % 'maxMinFrac threshold'
    NPhaseLabels(1:nP) = {legendarray('nPh = ', allNPhases )};           % 'maxMinFrac threshold'
    msPvalLabels(1:nP) = {legendarray('p < ', pvalLims(2:end) )};           % p-value
    
    colorCodingVarNames = {'F1/DC_pref pairing type', 'F1/DC_cmp pairing type', 'F1/DC threshold', 'maxMinFrac threshold', 'num phases', 'p-value',  '[none]'};
    colorCodingLegLabels = { F1oDCpairingLabels,        F1oDCpairingLabels,    F1oDCthresholdLabels, MinFracThresholdLabels, NPhaseLabels, msPvalLabels, {{' '}} };
    nCC_cats     = cellfun(@(C) length(C{1}), colorCodingLegLabels);
    
    colorCodings = struct('name', colorCodingVarNames, 'labels', colorCodingLegLabels, 'nCC_cats', nCC_cats);
    
    if ~useMeasurePvals
        colorCodings(find(strcmp(colorCodingVarNames, 'p-value'),1)) = [];
    end
    colorCodingVariable = colorCodings(1).name;  %'F1/DC pairing type',
    colorCodingVariable_prev = '';
    nCCs = length(colorCodings);
%     colorCodingNicknames = {'F1oDC_pairing', 'F1oDC_threshold', 'maxMinFrac_threshold', 'p_value', 'none'};

    % Legend variables   

    nMaxCats = max(nCC_cats(:))+1;  % add 1 so that colors look nicer
    colorMasks = cell(nLoc, nM, nMaxCats);
    
    % Color coding of different pair types
    h_ax = zeros(nP, nLoc, nM );
    h_bars = zeros(nP, nLoc, nM, nMaxCats );
    h_title = zeros(nP, nLoc, nM );
    h_figtitle = zeros(nP, 1);
    h_ylabel = zeros(nP, nLoc, nM );
    h_xlabel = zeros(nP, nLoc, nM );
    h_legend = zeros(nP);
    h_err    = zeros(nP, nLoc, nM); 
    [subM, subN, prevSubM, prevSubN] = deal([]);
    h_dummy_ax = zeros(nP,1);
    h_dummy_bars = zeros(nP, nMaxCats);
    h_legend = zeros(nP,1);

    h_sfPairLabel = zeros(nP, 1);    
    h_sfLocLabel = zeros(nLoc, 1);
    h_sfMsLabel = zeros(nLoc, nM);    
    h_sumFigBox = cell(nP, nLoc, nM);
    h_sfText = zeros(nP, nLoc, nM);
    h_sfLine = 0;
    h_sfFilterStr = 0;
    
    idxMtx = zeros(nUnits, nUnits, 'uint32');
    idxMtx(pairIdxList) = 1:length(pairIdxList);
%     idxMtx_ind = find(idxMtx(:));
    fontsize = 10;
    
    
    msVsMsFig = 28;
    [idx_ms_x, idx_ms_y, idx_ms_xy, idx_wcc, h_ms12_L, h_ms12_tit] = deal(0);
    [ms_x, ms_y] = deal('');
    
    sumFigId = 20+gratingType;
    sumFigSpacing = 2;
%     sumFigPairRows =  {'Wcc', 'Wrcc', 'Bcc',  'Wcm', 'Wrcm', 'Bcm'};
    sumFigAxes = 0;    

    firstTime = true;
    [curBin_N, curBin_dN, curSetVals, curSetIdx, curDistPval, curPvalLabels] = deal(cell(nP, nLoc, nM));    
    [curDistMean, curDistStd, curDistStderr, curDistMedian, curDistN] = deal(zeros(nP, nLoc, nM));        


    toDisplayNames = {'Histograms', 'Summary', 'MsVsMs'};
    toDisplayCur = [true, false, false];
    toDisplayPrev = toDisplayCur;    
    
    for j = 1:length(pairIdxs)
        if iscell(pairIdxs{j})
            pairIdxs{j} = cellfun(@(inds) idxMtx( inds ), pairIdxs{j}, 'un', false);
        else
            pairIdxs{j} = idxMtx( pairIdxs{j} );
        end
    end
    
        % calculate color coding masks once (so can access faster later on)                
        
    for cc_idx = 1:nCCs     
        nCat = nCC_cats(cc_idx);
        colorMasks = cell(nLoc, nM);
        switch colorCodings(cc_idx).name
            
            case 'F1/DC_pref pairing type', 
                colorMask = zeros(nPairs,1, 'uint8');
                SCtypes = [2, 0, 1.5, 1]; % [0, 1, 1.5, 2] = [cC-mC / cC-mS / cS-mC / cS-mS]
%                 SCtypes = [0, 1, 1.5, 2]; % [0, 1, 1.5, 2] = [cC-mC / cC-mS / cS-mC / cS-mS]
                for c_j = 1:nCat            % [0, 1,      2] = [CC,       SC/CS,          SS];
                    idxThisCategory = pairData.SCtype_pref == SCtypes(c_j);
                    colorMask(idxThisCategory) = c_j+10; 
                end
%                 assert(all(colorMask));
                colorMask = colorMask-10;
                colorMasks(:,:) = {colorMask};
                
            case 'F1/DC threshold', 
                colorMask = zeros(nPairs,1, 'uint8');
                for c_j = 1:nCat
                    boundStr = iff(c_j < nCat, '[)', '[]');
                    idxThisCategory = between(pairData.minF1oDC_pref, F1oDC_thresholds(c_j), F1oDC_thresholds(c_j+1), boundStr);
                    colorMask(idxThisCategory) = c_j;
                end
                assert(all(colorMask));
                colorMasks(:,:) = {colorMask};
                
            case 'maxMinFrac threshold',                
                loc_data = pairData.loc_minFracOfMaxes; 
                assert(all(size(loc_data) == [nPairs nLoc]));
                for loc_j = 1:nLoc
                    colorMask = zeros(nPairs,1, 'uint8');
                    for c_j = 1:nCat
                        boundStr = iff(c_j < nCat, '[)', '[]');                        
                        idxThisCategory = between(loc_data(:,loc_j), maxMinFrac_thresholds(c_j), maxMinFrac_thresholds(c_j+1), boundStr);
                        colorMask(idxThisCategory) = c_j;
                    end
%                     assert(all(colorMask | isnan(loc_data(:,loc_j))));
                    colorMasks(loc_j,:) = {colorMask};
                end

            case 'F1/DC_cmp pairing type', 
                loc_data = pairData.SCtype_cmp;
                SCtypes = [0, 1, 1.5, 2]; % [0, 1, 1.5, 2] = [cC-mC / cC-mS / cS-mC / cS-mS]
                for loc_j = 1:nLoc
                    colorMask = zeros(nPairs,1, 'uint8');
                    for c_j = 1:nCat            % [0, 1,      2] = [CC,       SC/CS,          SS];
                        idxThisCategory = loc_data(:,loc_j) == SCtypes(c_j);
                        colorMask(idxThisCategory) = c_j;
                    end
                    assert(all(colorMask | isnan(loc_data(:,loc_j))));
                    colorMasks(loc_j,:) = {colorMask};
                end
                
                
            case 'num phases',
                colorMask = zeros(nPairs,1, 'uint8');                
                for c_j = 1:nCat          % [0, 1,      2] = [CC,       SC/CS,          SS];
                    idxThisCategory = pairData.n_phases == allNPhases(c_j);
                    colorMask(idxThisCategory) = c_j;
                end
                assert(all(colorMask));
                colorMasks(:,:) = {colorMask};
                
            case 'p-value';
                for loc_j = 1:nLoc
                    for m_j = 1:nM                        
                        if isfield(S{loc_j, m_j}, 'pval')
                            colorMask = zeros(nPairs,1, 'uint8');
                            for c_j = 1:nCat                            
                                boundStr = iff(c_j < nCat, '[)', '[]');
                                idxThisCategory = between(S{loc_j, m_j}.pval, pvalLims(c_j), pvalLims(c_j+1), boundStr);
                                colorMasks{loc_j, m_j}(idxThisCategory) = c_j;
                            end
                        else
                            colorMask = ones(nPairs,1, 'uint8');
                        end
                        assert(all(colorMask));
                        colorMasks(loc_j, m_j) = {colorMask};
                    end
                end
                
            case '[none]',
                colorMasks(:,:) = {ones(nPairs,1, 'uint8')}; 
        end
        
        colorCodings(cc_idx).lists = colorMasks;        
    end


    [showHistograms, showSummaryFig, showMsVsMs] = dealV( toDisplayCur );
        
    histFirstTime = true;        
    summaryFirstTime = true;
    msVsMsFirstTime = true;
    
    plotOrUpdateAllFigures;    

    %%%%%%%%%%%%%%%%%%%%%  CONTROL PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if skipControlPanel
%         return;
%     end

    mainTitle = [titleCase(gratingType_s) ' Gratings'];
    panel_fig_id = 160 + gratingType*10;   
    
    
        
    function fprintf2(varargin)
        if printProgress, fprintf(varargin{:}); end
    end
        
    function s = cmpFunc2Str(cmpFunc)
        if isequal(cmpFunc, @ge),     s = '\geq';
        elseif isequal(cmpFunc, @le), s = '\leq';
        elseif isequal(cmpFunc, @lt), s = '<';
        elseif isequal(cmpFunc, @gt), s = '>';
        elseif isequal(cmpFunc, @eq), s = '=';
        end        
    end    
    
    function makeVisible(h)
        if iscell(h)
            cellfun(@makeVisible, h);
            return;
        end
        h = h(:);
        invisible = strcmp(get(h, 'visible'), 'off');
        if any(invisible)
            set(h(invisible), 'visible', 'on');
        end
    end


    
    function [distMean, distStd, distStderr, distN, distPval, bin_N, bin_dN, setVals, setIdx] = plotOrUpdateAllFigures
    % INITIALIZE - FIRST TIME PLOTTING
        [bin_N,    bin_dN,    distMean,    distStd,    distStderr,    distMedian,    distN,    distPval,   setVals,    setIdx] = deal(... 
         curBin_N, curBin_dN, curDistMean, curDistStd, curDistStderr, curDistMedian, curDistN, curDistPval, curSetVals, curSetIdx);
        pvalLabels = curPvalLabels;        
     
%         fprintf('Starting to recalculate: '); toc;        
        pairTypeInds = find(pairTypesActive);
        locationInds = find(locationsActive);
        measureInds = find(measuresActive);
        

%         nPairsEachType = cellfun(@length, allPairIdxs_list);
%         pairTypesActive_tmp =pairTypesActive; 
%         if pairTypesActive_tmp(3), pairTypesActive_tmp(2) = 0; end
%         if pairTypesActive_tmp(6), pairTypesActive_tmp(5) = 0; end
%         nPairsInUse = sum(nPairsEachType(pairTypesActive_tmp));
%         usePairFilter = nPairsInUse/nPairs < .2;
        
%         if usePairFilter
%             usePairFilter
%         end
        
        fprintf2('\nSTART:');
        [updateColorMasks, updateLegend] = deal(~strcmp(colorCodingVariable, colorCodingVariable_prev));
                        
        cc_i = find( strcmp( colorCodingVariable, {colorCodings.name}));
        useColorCoding = ~strcmp(colorCodingVariable, '[none]') || toDisplayCur(2);  % 2=summary

%         if useColorCoding            
            colorMasks = colorCodings(cc_i).lists;
%         end
        activePsLocMs = bsxfun(@and, bsxfun(@and, pairTypesActive', locationsActive), reshape(measuresActive, [1 1 nM]));        
        activeLocMs = bsxfun(@and, locationsActive', measuresActive);                    

        
        %%%%%% 1. UPDATE DATA FILTERS & MASKS               
        % update if: (1) changed from active/inactive or (2) is active, and threshold value has changed.         

        to4D = @(x) reshape(x, [1 1 1 length(x)]);        
        filterActive4 = to4D(filterActive);
        setFilterActive_prev4 = cell2mat(cellfun(to4D, setFilterActive_prev, 'un', 0));
                
%         th_to4D = @(s_th) to4D( cellfun(@(fv1) double(s_th.(fv1)), filterFields));
        th_to4D = @(s_th) to4D( cell2mat(struct2cell(s_th)) );                
        th4 = th_to4D(th); 
        setTh_prev4 = cell2mat(cellfun(th_to4D, setTh_prev, 'un', 0));
                        
        setsWithFiltersChanged = bsxfun(@and, activePsLocMs,  bsxfun(@ne, filterActive4, setFilterActive_prev4));        
        setsWithFiltersActiveAndDiffTh = bsxfun(@and, filterActive4, bsxfun(@ne, th4, setTh_prev4));
        setsWithFiltersToUpdate = bsxfun(@and, activePsLocMs, (setsWithFiltersChanged | setsWithFiltersActiveAndDiffTh));
        setsToUpdate = (any(setsWithFiltersToUpdate, 4));
        locMsToUpdate = reshape(any(setsToUpdate, 1), [nLoc, nM]) | firstTime;
        
        
        
        % 1. update filters that pertain only on the pairs of cells, not on location / measure.
        pairThToUpdate = squeeze(any(any(any(setsWithFiltersToUpdate(:,:,:,pairFilterInds), 1),2),3))';        
        activePairThToUpdate = pairThToUpdate & filterActive(pairFilterInds);

        pairsWithFiltersChanged = setsWithFiltersChanged(:,:,:,pairFilterInds);  % (note: dont want to use pairsWithFiltersToUpdate here (which includes if FiltersActiveAndDiffTh): don't care if newly added loc/measure wasn't present: can still use same filters.)
        updatePairMask = any(pairsWithFiltersChanged(:)); % also update pair mask if a pair filter was turned on/off, even if no thresholds were changed
        
        if (nnz(activePairThToUpdate > 0))             
            for pf_i = find(pairThToUpdate & filterActive(pairFilterInds)) % (a) update pair masks that are active and have changed 
                fltr = pairFilters(pf_i);
                filterVar = fltr.filterName; th_name = [filterVar '_th'];
                if ~isfield(pairSubMasks, th_name) || pairSubMasks.(th_name) ~= th.(filterVar);
                    pairSubMasks.([filterVar '_mask']) = fltr.op(pairData.(fltr.varName), th.(filterVar) ); 
                    pairSubMasks.(th_name) = th.(filterVar);
                    fprintf2('P[%s]',filterVar(1:5));
                    updatePairMask = true;
                end
            end
        end            
        
        if updatePairMask
            pairMask = true(nPairs, 1);    % (b) merge all pair masks that are active (if one has been changed)
            for pf_i = find(filterActive(pairFilterInds))
                filterVar = pairFilters(pf_i).filterName;
                pairMask = pairMask & pairSubMasks.([filterVar '_mask']);
                fprintf2('P-M[%s]',filterVar(1:5));                                            
            end            
        end

        % 2. Update filters that are depend on location        
        locThToUpdate = reshape(any(any(setsWithFiltersToUpdate(:,:,:,locFilterInds), 1),3), [nLoc length(locFilterInds)]);        
        activeLocThToUpdate = bsxfun(@and, locThToUpdate, filterActive(locFilterInds) );

        locsWithFiltersChanged = setsWithFiltersChanged(:,:,:,locFilterInds);
        updateLocMasks = locationsActive.*any(locsWithFiltersChanged(:)); % also update all (active) location mask if a pair filter was turned on/off, even if no thresholds were changed;

        
        if (nnz(activeLocThToUpdate > 0))                       
            for loc = find(any(locThToUpdate,2)')        % (a) update location masks that have changed and need updating
                for lf_i = find(locThToUpdate(loc,:) & filterActive(locFilterInds))  
                    fltr = locFilters(lf_i);
                    filterVar = fltr.filterName; th_name = [filterVar '_th'];
                    if ~isfield(locSubMasks{loc}, th_name) || locSubMasks{loc}.(th_name) ~= th.(filterVar);
                        loc_data = pairData.(['loc_' fltr.varName]);
                        locSubMasks{loc}.([filterVar '_mask']) = fltr.op(loc_data(:,loc), th.(filterVar) );
                        locSubMasks{loc}.(th_name) = th.(filterVar);
                        fprintf2('L[%s:%s]',locationLabels{loc}(end-4:end),filterVar(1:5));
                        updateLocMasks(loc) = true;
                    end
                end
            end
        end
        for loc = find(updateLocMasks)   % (b) merge all location masks that are active and have had at least one of their filters updated
            locMasks{loc} = true(nPairs, 1);
            for lf_i = find(filterActive(locFilterInds))
                filterVar = locFilters(lf_i).filterName;
                locMasks{loc} = locMasks{loc} & locSubMasks{loc}.([filterVar '_mask']);
                fprintf2('L-M[%s:%s]',locationLabels{loc}(end-4:end),filterVar(1:5));
            end            
        end
        
        % 3. Update filters that depend on location and measure
%         valThToUpdate = squeeze(any(setsWithFiltersToUpdate(:,:,:,[10 11]), 1));
        valThToUpdate = reshape(any(setsWithFiltersToUpdate(:,:,:,valFilterInds), 1), [nLoc, nM]);
        activeValThToUpdate = bsxfun(@and, valThToUpdate, filterActive(valFilterInds) );

        valsWithFiltersChanged = setsWithFiltersChanged(:,:,:,valFilterInds);
        updateValMasks = activeLocMs.*any(valsWithFiltersChanged(:)); % also update all (active) val masks if a pair filter was turned on/off, even if no thresholds were changed;
        
        if (nnz(activeValThToUpdate) > 0)
            for loc = find(any(any(valThToUpdate,2),3)')                                
                for ms = find(any(valThToUpdate(loc,:,:),3))                    
                    for vf_i = find(valThToUpdate(loc,ms,:) & filterActive(valFilterInds))'
                        fltr = valFilters(vf_i);
                        filterVar = fltr.filterName; th_name = [filterVar '_th'];
                        if ~isfield(valSubMasks{loc, ms}, th_name) || valSubMasks{loc, ms}.(th_name) ~= th.(filterVar);
                            if strcmp(filterVar, 'pval') && ~isfield(S{loc, ms}, 'pval'), continue, end;                        
                            valSubMasks{loc, ms}.([filterVar '_mask']) = fltr.op([S{loc, ms}.(fltr.varName)], th.(filterVar) );
                            valSubMasks{loc, ms}.(th_name) = th.(filterVar);
%                             fprintf2('V[%s:%s:%s]',locationLabels{loc}(end-4:end),measureLabels{ms}(end-1:end),filterVar(1:4));                            
                            updateValMasks(loc,ms) = true;
                        end
                    end
                end
            end            
        end
        for loc = find(any(updateValMasks,2)')
            for ms = find(updateValMasks(loc,:))
                valMasks{loc, ms} = true(nPairs, 1);
                for vf_i = find(filterActive(valFilterInds))
                    filterVar = valFilters(vf_i).filterName;
                    if strcmp(filterVar, 'pval') && ~isfield(S{loc, ms}, 'pval'), continue, end;
                    valMasks{loc, ms} = valMasks{loc, ms} & valSubMasks{loc, ms}.([filterVar '_mask']);
%                     fprintf2('V-M[%s:%s:%s]',locationLabels{loc}(end-4:end),measureLabels{ms}(end-1:end),filterVar(1:4));
                end
            end
        end

        
        for loc = find(any(locMsToUpdate,2))' 
            for ms = find(locMsToUpdate(loc,:))
                fprintf2('Merge[%s:%s]',locationLabels{loc}(end-4:end),measureLabels{ms}(end));
                dataMasks{loc, ms} = pairMask & locMasks{loc} & valMasks{loc, ms};
            end
        end
        
        
%         filterStr = filterActive
        
        %%%%%% 2. Merge filters, calculate occupancy of each bin
         % for errorbars (for resampled bins)
                
        nSummedPhis = switchh(opt.defaultPhase_oeState, {'aa', 'oe'}, [1, 2]);
        nResamplesChange = false(nP,1); 
        nResamplesChange(cellfun(@(x) ~isempty(x), strfind(pairTypes, '_r'))) = nResamples ~= nResamplesPrev; % only relevant for pairs with resampling.
%         displayChange = ~strcmp(toDisplay, toDisplayPrev) && ~strcmp(toDisplayPrev, 'Both');
        displayChange = any(~toDisplayPrev & toDisplayCur);
        
        psLocMsGraphsToUpdate = bsxfun(@or, setsToUpdate | firstTime | displayChange | (updateColorMasks && showHistograms ), nResamplesChange); 
        locMsGraphsToUpdate = reshape(any(psLocMsGraphsToUpdate,1), [nLoc, nM]);               
        
        nCat = nCC_cats(cc_i);
  
        rescaleDphiDistribution = 1;
        rescaleBccDistribution = 1;
        
        fprintf2('\n');
        for loc = find(any(locMsGraphsToUpdate,2))'
            for ms = find(locMsGraphsToUpdate(loc,:))

                is_pval = isfield(S{loc, ms}, 'pval');
                is_dphi = any(strcmp(measures{ms}, {'dphi', 'dF1'})) && ~strncmp(locations{loc}, 'wgtedSum', 7);
                
%                 binMask = S{loc, ms}.binMask;
                binIds = S{loc, ms}.bins;
                
                binEdges = S{loc, ms}.binEdges;
                nBins = length(binEdges)-1;
                                
                useColorCodingHere = useColorCoding;
                if strcmp(colorCodingVariable , 'p-value') && ~is_pval % some measures don't have p-values.
                    useColorCodingHere = false;
                end
                                                
                % Calculate occupancy of each bin.                
                mask = dataMasks{loc, ms};
                                
%                 if useColorCodingHere
                    catIds = colorMasks{loc, ms};
%                 else
%                     catIds = [];
%                 end
                nph_idx = find(strcmp('num phases', {colorCodings.name}),1);                
                nPhaseMask = colorCodings(nph_idx).lists{1};                 
                                                
                for pr = find(squeeze(psLocMsGraphsToUpdate(:,loc,ms))')
                    
                   
                    fprintf2('Up[%s:%s:%s]',locationLabels{loc}(end-4:end),measureLabels{ms}(end), pairTypes{pr});
                                        
                    pairVecIdxs = pairIdxs{pr};
                    useContrib = false;%strncmp(pairTypes{pr}, 'Wcc', 3);
                    if ~iscell(pairVecIdxs)  % ie. not randomized pairs.
%                         fprintf2('%d %d %d ', loc, ms, pr)                        

%                         rescaleDphi_now = rescaleDphiDistribution && is_dphi
%                         p_wcc = find(strcmp(pairTypes, 'Wcc'), 1);  nTotalWccBins = sum(bin_N{p_wcc, loc, ms}(:));
                        rescaleBcc_now = rescaleBccDistribution && any( strcmp('Wcc', pairTypes(pairTypesActive) )) && any( strcmp('Bcc', pairTypes(pairTypesActive) ));
                                                                                                    
                        
                        if rescaleBcc_now                            
                            catIds_here = nMaxCats*(nPhaseMask-1);  % color coding categories are rescaled, as well as # phases, 
                            if ~isempty(catIds)                     % so temporarily pretend they are all different categories
                                catIds_here = catIds_here + catIds;
                            end
                            nMaxCats_here = length(allNPhases)*nMaxCats;
                        else 
                            catIds_here = catIds;
                            nMaxCats_here = nMaxCats;
                        end

%                         [bin_N{pr, loc, ms}, setVals{pr,loc,ms}, setIdx{pr,loc,ms}] = binCountForPairs_Matlab(pairIdxs{pr}, mask, binIds, catIds_here, nBins, nMaxCats_here, S{loc, ms}.val);                                                  
                        [bin_N{pr, loc, ms}, setVals{pr,loc,ms}, setIdx{pr,loc,ms}] = binCountForPairs_c(     pairIdxs{pr}, mask, binIds, catIds_here, nBins, nMaxCats_here, S{loc, ms}.val);
    
                        assert( isequal( S{loc, ms}.val( setIdx{pr,loc,ms}), setVals{pr,loc,ms}' ))
                        %                          bin_N_beforeRescaling{pr, loc, ms} = bin_N{pr, loc, ms};
                         
                        
                        if all([pr, loc, ms] == [2 1 1])
                            3;
                        end

                    
                        
                        if rescaleBcc_now && strcmp(pairTypes{pr}, 'Bcc') 
                            if rescaleDphiDistribution && is_dphi
                                3;
                            end
                            bcc_bins = bin_N{pr, loc, ms};
                            bcc_sumBins = sum(bcc_bins,1);
                            nTotalBccBins = sum(bcc_sumBins(:));
                            
                            p_wcc = find(strcmp(pairTypes, 'Wcc'));                                                        
                            wcc_bins = bin_N{p_wcc, loc, ms};
                            wcc_sumBins = sum(wcc_bins,1);
                            nTotalWccBins = sum(wcc_sumBins);
                            
                            if nTotalWccBins > 0   % if nothing in wcc bins, can't use them to rescale
                                binRescale = zeros(size(wcc_sumBins));
                                b_idx = (wcc_sumBins > 0 & bcc_sumBins);
                                binRescale(b_idx) = (wcc_sumBins(b_idx)./ bcc_sumBins(b_idx)) * (nTotalBccBins / nTotalWccBins);

                                for i = find(b_idx)
                                    bcc_bins(:,i) = bcc_bins(:,i) .* binRescale(i);
                                end
                                bcc_bins = bcc_bins * nTotalBccBins / sum(bcc_bins(:)); % rescale back to original numbers.

                                if rescaleDphiDistribution && is_dphi
                                    [setNPhases_wcc, setNPhaseCounts_wcc] = uniqueCount([pairData.n_phases(setIdx{p_wcc,loc,ms})]);
                                    [dPhiNull_val_wcc, dPhiNull_Count_wcc] = deltaPhiNull(setNPhases_wcc, setNPhaseCounts_wcc, allNPhases);
                                    dPhiNull_count_wcc_binned = rebinnedDPhiCount(dPhiNull_val_wcc, dPhiNull_Count_wcc, dPhi_binE);

                                    [setNPhases_bcc, setNPhaseCounts_bcc] = uniqueCount([pairData.n_phases(setIdx{pr,loc,ms})]);
                                    [dPhiNull_val_bcc, dPhiNull_Count_bcc] = deltaPhiNull(setNPhases_bcc, setNPhaseCounts_bcc, allNPhases);
    %                                 dPhiNull_count_bcc_binned = rebinnedDPhiCount(dPhiNull_val_bcc, dPhiNull_Count_bcc, dPhi_binE);

                                    % we have already scaled the bcc distribution to look like the wcc distribution, so divide by the 
                                    % null number acc to wcc npairs ratio (but scaled by ratio of n(bcc)/n(wcc) 
                                    n_ratio = sum(dPhiNull_Count_bcc)/sum(dPhiNull_Count_wcc);
                                    bcc_bins = bsxfun(@rdivide, bcc_bins,  dPhiNull_count_wcc_binned(:)*n_ratio);
                                    wcc_bins = bsxfun(@rdivide, wcc_bins,  dPhiNull_count_wcc_binned(:));
                                end
                            end
                            
                            wcc_bins = reshape(wcc_bins, [nBins, nMaxCats, length(allNPhases)]);
                            bin_N{p_wcc, loc, ms} = sum(wcc_bins,3);
                            bcc_bins = reshape(bcc_bins, [nBins, nMaxCats, length(allNPhases)]);
                            bin_N{pr, loc, ms} = sum(bcc_bins,3);
    %                             if any(bcc_sumBins > 0)
    %                                 assert( abs( sum(bcc_bins(:))-nTotalBccBins) < 1e-10)                            
    %                             end
                            
    
                        elseif ~rescaleBcc_now && rescaleDphiDistribution && is_dphi                            
                            [setNPhases, setNPhaseCounts] = uniqueCount([pairData.n_phases(setIdx{pr,loc,ms})]);
                            [dPhiNull_val, dPhiNull_Count] = deltaPhiNull(setNPhases, setNPhaseCounts, allNPhases, nSummedPhis);
                            dPhiNull_count_binned = rebinnedDPhiCount(dPhiNull_val, dPhiNull_Count, dPhi_binE);
                            binCounts = bin_N{pr, loc, ms};                            
                            rescaledBinCounts = bsxfun(@rdivide, binCounts,  dPhiNull_count_binned(:) +eps);
                            bin_N{pr, loc, ms} = rescaledBinCounts;
                        end
                        if all(isnan(bin_N{pr, loc, ms}(:)))
                            3;
                        end
                           
                            
                        takeUniqueCell1 = false;
                        if strcmp(pairTypes{pr}, 'Bcc') && takeUniqueCell1
                            
                            Gs = cat(1, pairData.Gids(setIdx{pr,loc,ms}) );
                            Cs = cat(1, pairData.cellIds(setIdx{pr,loc,ms}) );
                            
                            cellIdxs1 = arrayfun(@(G, C) find([allCells.Gid] == G & [allCells.cellId] == C), Gs(:,1), Cs(:,1));
                            cellIdxs2 = arrayfun(@(G, C) find([allCells.Gid] == G & [allCells.cellId] == C), Gs(:,2), Cs(:,2)); 
                            
                            OIidx = cat(2, pairData.oriSp_cmp(setIdx{pr,loc,ms}))';
                            GCOI1 = [Gs(:,1), Cs(:,1), OIidx];
                            GCOI2 = [Gs(:,2), Cs(:,2), OIidx];
%                             GCOI = [Gs(:,1), Cs(:,1), OIidx(:,2)];
                            [uGCOI1, GCOI_list1] = uniqueList(GCOI1, 'rows');
                            GCOI_count1 = cellfun(@length, GCOI_list1);
                            GCOI_list1_1st = cellfun(@(x) x(1), GCOI_list1);
                            [uGCOI2, GCOI_list2] = uniqueList(GCOI2, 'rows');
                            GCOI_count2 = cellfun(@length, GCOI_list2);
                            GCOI_list2_1st = cellfun(@(x) x(1), GCOI_list2);
                            
                            oi_select = zeros(36,10);
                            for i = 1:length(GCOI1)
                                oi_select(OIidx(i,1), OIidx(i,2)) = oi_select(OIidx(i,1), OIidx(i,2)) + 1;
                            end
                            figure(3); imagesc(oi_select); colorbar;
                            title({'ori/sp selected ', sprintf('(N = %d)', length(GCOI1))});
                            title({'ori/sp selected (#2 unique) ', sprintf('(N = %d)', length(GCOI_list2_1st) )})
                            

                            tc = [pairData.phaseTCs(setIdx{pr,loc,ms})];
                            tc1 = tc(:,:,1);
                            tc2 = tc(:,:,2);
                            for i = 1:size(tc1, 2);
                                tc1(:,i) = tc1(:,i) / norm(tc1(:,i));
                                tc2(:,i) = tc2(:,i) / norm(tc2(:,i));
                            end
                            
                            meanTc1 = mean(tc1, 2);
                            meanTc2 = mean(tc2, 2);
                            meanUTc1 = mean(tc1(:,GCOI_list2_1st), 2);
                            meanUTc2 = mean(tc2(:,GCOI_list2_1st), 2);

                            [aa, bb] = uniqueList( tc1(:,GCOI_list1_1st)', 'rows');
                            
                            
                            figure(55)
                            ph = linspace(0, 360, length(meanTc1)+1); ph = ph(1:end-1);
                            plot(ph, meanTc1, 'bo-', ph, meanTc2, 'gs:');
                            set(gca,'xtick', ph)
                            3;
                            plot(ph, meanUTc1, 'bo-', ph, meanUTc2, 'gs:');
                            
                            
%                             setIdx{pr,loc,ms} = GCOI_list1;
%                             setVals{pr,loc,ms} = S{loc,ms}.val( setIdx{pr,loc,ms});                                                       
                            
                        end
                        
                        
                        
%                         fprintf('ok\n');
                        if useContrib                            
                            distMean(pr, loc, ms) = wgt_nanmean(setVals{pr,loc,ms}, curcontribs);
                            distStd(pr, loc, ms) =  wgt_nanstd(setVals{pr,loc,ms}, curcontribs);
                            distStderr(pr, loc, ms) =  wgt_nanstderr(setVals{pr,loc,ms}, curcontribs);                            
                        else
                            distMean(pr, loc, ms) = nanmean(setVals{pr,loc,ms}); 
                            distStd(pr, loc, ms) =  nanstd(setVals{pr,loc,ms});
                            distStderr(pr, loc, ms) =  nanstderr(setVals{pr,loc,ms});
                            distMedian(pr, loc, ms) =  nanmedian(setVals{pr,loc,ms});
                        end                  
                        distN(pr, loc, ms) = nnz(~isnan(setVals{pr,loc,ms}));
                        
                        
                        
                    else
                        Ni = zeros(nBins, nMaxCats, nResamples);
                        sampleVals = cell(1, nResamples);
%                         fprintf('%d %d %d ', loc, ms, pr)
                        for i = 1:nResamples
%                             fprintf('(%d)', i)
%                             [Ni(:,:,i), sampleVals{i}] = binCountForPairs_Matlab(pairVecIdxs{i}, mask, binIds, catIds, nBins, nMaxCats, S{loc, ms}.val);
                            [Ni(:,:,i), sampleVals{i}] = binCountForPairs_c(pairIdxs{pr}{i}, mask, binIds, catIds, nBins, nMaxCats, S{loc, ms}.val);
                        end
%                         fprintf('ok\n');
                        bin_N{pr, loc, ms} = mean(Ni, 3);                        
                        bin_dN{pr, loc, ms} = std(sum(Ni,2), 0, 3); % sum across categories, then take std across samples.                                                     
                        
                        sampleMeans = cellfun(@nanmean, sampleVals);
                        distMean(pr, loc, ms)   = mean(sampleMeans); 
                        distStd(pr, loc, ms)    = mean(cellfun(@nanstd, sampleVals)); %, contrib));
                        distStderr(pr, loc, ms) = mean(cellfun(@nanstderr, sampleVals)); %, contrib));
                        distN(pr, loc, ms)      = mean(cellfun(@(v) nnz(~isnan(v)), sampleVals));
                    end                    
                    
                    if any(isnan(bin_N{pr, loc, ms}(:)))
                        3;
                    end
                end
                3;
%                 nMaxCats
%                 assert( size(bin_N
                
                
            end
        end
        
%         xallNs = [bin_N{:}];
%         if any(isnan(xallNs(:)))
%             3;
%         end

        %%%%%% 3. Perform significance tests, calculate (p-values)
        p_idx = [];
        doTtests = true;  
            do1SmpTtestAlways = 1;
        doMedianTests = false;
        doUtests = true;
        doKStests = false;
        doChiSqrTests = true;
            doChiSqrSeparate = true;
            doChiSqrCombined = false;
        doF1F2tests = false;
        
        for ms = find(any(locMsGraphsToUpdate,1)) 
            
            isCC = ~isempty(strfind(measures{ms}, 'cc'));
            for loc = find(locMsGraphsToUpdate(:,ms))'
                for pr = find(squeeze(psLocMsGraphsToUpdate(:,loc,ms)))'
                                        
                    p_idx = 1;
                    distPval{pr, loc, ms} = nan(1, nMaxSigTests);
                    pvalLabels{pr,loc,ms} = {};
                    
                    if isempty(setVals{pr,loc,ms})
                        continue;
                    end                    
                    
                    if any( strcmp(pairTypes{pr}, {'Wcc', 'Wcm'}) )
                        pairTypeCmp = strrep(pairTypes{pr}, 'W', 'B');
                        p_cmp = find(strcmp(pairTypeCmp, pairTypes));
                        if (~isempty(p_cmp) && isempty(setVals{p_cmp,loc,ms})), p_cmp = []; end;
                    else
                        p_cmp = [];
                    end
                    
                    [setNPhases, setNPhaseLists] = uniqueList(pairData.n_phases(setIdx{pr,loc,ms}));
                    setNPhaseCounts = cellfun(@length, setNPhaseLists);
                    set_nph = {setNPhases, setNPhaseLists};
                    if ~isempty(p_cmp)
                        [setNPhases_cmp, setNPhaseLists_cmp] = uniqueList(pairData.n_phases(setIdx{p_cmp,loc,ms}));
                        set_nph_cmp = {setNPhases_cmp, setNPhaseLists_cmp};                        
                    end
                    
                    % 1. MW-Wilcoxon-test (U-test).
                    if doUtests && ~isempty(p_cmp) && any(strcmp(measures{ms}, {'dot', 'cc', 'rho', 'STA_cc', 'MID_cc'})) % all except dphi
                                                
                        if strcmp(measures{ms}, 'dF1')  % don't separate phases
                            distPval{pr, loc, ms}(p_idx) = myRanksum(setVals{pr,loc,ms}, setVals{p_cmp,loc,ms});  % 'right'
                            pvalLabels{pr,loc,ms}{p_idx} = 'p_{U}'; p_idx = p_idx+1; 
                        else                            % separate phases
                            distPval{pr, loc, ms}(p_idx) = sectionedSigTest(...
                                'myRanksum', setVals{pr,loc,ms}, set_nph, setVals{p_cmp,loc,ms}, set_nph_cmp); % 'right'
                            pvalLabels{pr,loc,ms}{p_idx} = 'p_{U}'; p_idx = p_idx+1;                                                        
                        end
%                         distPval{pr, loc, ms}(p_idx) = myRanksum(setVals{pr,loc,ms}, setVals{p_cmp,loc,ms});  % 'right'
%                         pvalLabels{pr,loc,ms}{p_idx} = 'p_{Uc}'; p_idx = p_idx+1; 
                        
                    end
                    
                    
                    % 2. T-test.
                    if doTtests % || isempty(p_cmp)
                        if any(strcmp(measures{ms}, {'cc', 'rho', 'dphi', 'dF1'})) || isCC
                            if isempty(p_cmp) || do1SmpTtestAlways                                
                                [~, distPval{pr, loc, ms}(p_idx)] = ttest(setVals{pr,loc,ms}, mMullMeans(ms)); 
                                pvalLabels{pr,loc,ms}{p_idx} = ['p_t'];  p_idx = p_idx+1;
                            end                            
                            if ~isempty(p_cmp)                            
                                [~, distPval{pr, loc, ms}(p_idx)] = ttest2(setVals{pr,loc,ms}, setVals{p_cmp,loc,ms}); 
                                pvalLabels{pr,loc,ms}{p_idx} = ['p_{t2}'];  p_idx = p_idx+1;
                            end
                        end
                    end
                                        
                    % 3. Wilcoxon (median) test.
                    if doMedianTests && any(strcmp(measures{ms}, {'cc', 'rho', 'dphi', 'dF1'}))  % this is not, in fact a good test for rho - ties make the null distribution asymmetric.
                        distPval{pr, loc, ms}(p_idx) = signrank(setVals{pr,loc,ms}, mMullMeans(ms));
                        pvalLabels{pr,loc,ms}{p_idx} = ['p_W']; p_idx = p_idx+1;
                    end
                    
                        
                    % 4. KS test.
                    if doKStests && any(strcmp(measures{ms}, {'dot', 'cc', 'rho',  'dF1'})) % all except dphi
                        
                        switch measures{ms}
                            case {'dot', 'cc', 'rho'}
                                 if ~isempty(p_cmp)   % compare with Bcc (after dividing phases)
                                    distPval{pr, loc, ms}(p_idx) = sectionedSigTest(...
                                        'kstest2', setVals{pr,loc,ms}, set_nph, setVals{p_cmp,loc,ms}, set_nph_cmp); % .01, 'larger'?
                                    pvalLabels{pr,loc,ms}{p_idx} = 'p_{ks}'; p_idx = p_idx+1;
                                 end
                                 
                            case 'dF1'
                                if ~isempty(p_cmp)  % compare with Bcc (without dividing phases)
                                    [tmp, distPval{pr, loc, ms}(end+1)] = kstest2(setVals{pr,loc,ms}, setVals{p_cmp,loc,ms});                                    
                                else                % compare with flat distribution
                                    dF1_CDF = [linspace(0, 180, 5)', linspace(0,1,5)'];
                                    [tmp, distPval{pr, loc, ms}(p_idx)] = kstest(setVals{pr,loc,ms}, dF1_CDF);                                    
                                end
                                pvalLabels{pr,loc,ms}{p_idx} = 'p_{ks}'; p_idx = p_idx+1;                                
                        end
                    end

                            
                    % 5. Chi-squared ) test.                                                        
                    if doChiSqrTests && any(strcmp(measures{ms}, {'dphi', 'dF1'})) && strcmp(cmpType_s, 'phase') ...
                            && ~strncmp(locations{loc}, 'wgtedSum', 7)
                        
%                             any( strcmp(pairTypes{pr}, {'Wcc', 'Wcm'}) )
%                             && ~isempty(p_cmp)
%                         if strcmp(measures{ms}, 'dphi') 
                            setNPhases = setNPhases(:)';
                            setNPhaseCounts = setNPhaseCounts(:)';
                            if doChiSqrCombined
                                [dPhiNull_val, dPhiNull_Count] = deltaPhiNull(setNPhases, setNPhaseCounts, [], nSummedPhis);
                                distPval{pr, loc, ms}(p_idx) = histChiSqrTest(setVals{pr,loc,ms}, dPhiNull_val, dPhiNull_Count);
                                pvalLabels{pr,loc,ms}{p_idx} = 'p_{\chic}';  p_idx = p_idx+1;
                            end
                            if doChiSqrSeparate 
                                [dPhiNull_vals, dPhiNull_Counts] = arrayfun(@(nph, cnt) deltaPhiNull(nph, cnt, [], nSummedPhis), setNPhases, setNPhaseCounts, 'un', 0);
                                nullDists = [num2cell(setNPhases); dPhiNull_vals; dPhiNull_Counts];
                                distPval{pr, loc, ms}(p_idx) = sectionedSigTest('chiSqr', setVals{pr,loc,ms}, set_nph, nullDists);
                                pvalLabels{pr,loc,ms}{p_idx} = 'p_{X}';  p_idx = p_idx+1;

                            end
%                         elseif 00 && strcmp(measures{ms}, 'dF1')                             
%                             dF1_binIds = S{loc,ms}.bins(setIdx{pr,loc,ms});                            
%                             nF1_bins = length(dF1_binE)-1;
%                             dF1null_vals = 1:nF1_bins;
%                             dF1null_N    = ones(1,nF1_bins) * length(dF1_binIds)/nF1_bins;
%                             distPval{pr, loc, ms}(p_idx) = histChiSqrTest(dF1_binIds, dF1null_vals, dF1null_N);
%                             pvalLabels{pr,loc,ms}{p_idx} = 'p_X'; p_idx = p_idx+1;
%                         end
                    end
                            
                                                            
                    if doF1F2tests && any(strcmp(measures{ms}, {'dphi', 'dF1'}))
                        if strcmp(measures{ms}, 'dphi')
                            binC = binEdge2cent(dPhi_binE);
                        elseif strcmp(measures{ms}, 'dF1')
                            binC = binEdge2cent(dF1_binE);
                        end

                        distPval{pr, loc, ms}(p_idx) = testSignificanceOfHistF1oDC(binC, sum(bin_N{pr,loc,ms},2), distN(pr,loc,ms), 1);
                        pvalLabels{pr,loc,ms}{p_idx} = 'p_{F1}'; p_idx = p_idx+1;
                        distPval{pr, loc, ms}(p_idx) = testSignificanceOfHistF1oDC(binC, sum(bin_N{pr,loc,ms},2), distN(pr,loc,ms), 2);
                        pvalLabels{pr,loc,ms}{p_idx} = 'p_{F2}'; p_idx = p_idx+1;
                        
                        % if ~isempty(p_cmp) -- should ideally compare with Bcc, not uniform
%                             [dPhiNull_vals, dPhiNull_Counts] = arrayfun(@(nph, cnt) deltaPhiNull(nph, cnt), setNPhases, setNPhaseCounts, 'un', 0);
%                             nullDists = [num2cell(setNPhases); dPhiNull_vals; dPhiNull_Counts];                            
%                                 [dPhiNull_vals_all, dPhiNull_Count_all] = deltaPhiNull(setNPhases, setNPhaseCounts, allNPhases);
%                                 dPhiNull_count_binned = rebinnedDPhiCount(dPhiNull_vals_all, dPhiNull_Count_all, dPhi_binE);                                

                    end
                            
                    
                end
            end
        end
        tmp1 = p_idx; %#ok<NASGU>

        showSigTestsOnPlots = true;
        showNOnPlots = true;
        fontsize = 9;
        
%         fprintf('Finished recalculating: '); toc;
                

        % PLOT results
        if showControlPanel && showHistograms
            skipYlabels = 0;
            gapAbove = iff(skipBigHeading, 0, .05);
            gapOnRight = iff(hideLegend, 0, .1);
            M_spacing = [gapAbove .01 .04];
            N_spacing = [0.1     .01 gapOnRight];  
            subM = length(measureInds);   % 4 2 subplotGap
            subN = length(locationInds);
            
            adjustAxesLocation = ~isequal([subM subN], [prevSubM prevSubN]);            
            
            % show/hide active/inactive figures
            isVisible = @(h) strcmp(get(h, 'visible'), 'on');
            axesVisible = arrayfun(isVisible, h_ax);
            axesToHide = h_ax & axesVisible & ~activePsLocMs;
            axesToShow = h_ax & ~axesVisible & activePsLocMs;

            barsVisible = arrayfun(isVisible, h_bars);
            barsToHide = bsxfun(@and, h_bars & barsVisible, ~activePsLocMs);
            barsToShow = bsxfun(@and, h_bars & ~barsVisible, activePsLocMs);

            ylabelsVisible = arrayfun(isVisible, h_ylabel);
            ylabelsNeed = false(nP, nLoc, nM);  ylabelsNeed( :, find(locationsActive, 1, 'first') , measuresActive ) = true;                
            ylabelsToHide = h_ylabel & ylabelsVisible & ~ylabelsNeed;
            ylabelsToShow = h_ylabel & ~ylabelsVisible & ylabelsNeed;

            xlabelsVisible = arrayfun(isVisible, h_xlabel);
            xlabelsNeed = false(nP, nLoc, nM);	xlabelsNeed( :, locationsActive, find(measuresActive, 1, 'last') ) = true;
            xlabelsToHide = h_xlabel & xlabelsVisible & ~xlabelsNeed;
            xlabelsToShow = h_xlabel & ~xlabelsVisible & xlabelsNeed;
            
            cellfun(@(hs) set(hs, 'visible', 'on'), ...
                {h_ax(axesToShow), h_bars(barsToShow), h_ylabel(ylabelsToShow), h_xlabel(xlabelsToShow), ...
                h_title(activePsLocMs), h_err(activePsLocMs & h_err) });
            cellfun(@(hs) set(hs, 'visible', 'off'), ...
                {h_ax(axesToHide), h_bars(barsToHide), h_ylabel(ylabelsToHide), h_xlabel(xlabelsToHide), ...
                h_title(~activePsLocMs), h_err(~activePsLocMs & h_err) });            
            
            for pr = pairTypeInds
                if histFirstTime
                    figure(pr+gratingType*10);
                    set(gcf, 'Name', [gratingType_s ' : ' pairTypes{pr}]);
                    clf;
                end

                % figure out best limits - make so that all measures have same y-axes limits             
                m_ylim2 = zeros(nM,2);
                loc_maxes = cell(1, length(locationInds));
                for ms = measureInds
                    [loc_maxes{:}] = bin_N{pr, locationInds, ms};
                    mx_data = max( cellfun(@(n) max(sum(n,2)), loc_maxes) );                    
                    m_ylim2(ms,:) = bestLimsFromData([0 mx_data]);
                end
                
                for loc_i = 1:length(locationInds)            
                    loc = locationInds(loc_i);
                    for m_i = 1:length(measureInds)
                        ms = measureInds(m_i);
    %                     activeCats = activeCats | any(bin_N,1);                        
                        pval_str = makePvalStr(distPval{pr, loc, ms}, pvalLabels{pr,loc,ms}); 
                        %pval_str = iff(any(strcmp(pairTypes{pr}, {'Wcc', 'Bcc'} )), ['. ' str_pval(distPval(pr, loc, ms))], '');
%                         showNOnPlots = false;
                        N_ifShow = iff(showNOnPlots, distN(pr, loc, ms), []);
%                         mnStd_str = str_mn_std(distMean(pr, loc, ms), distStderr(pr, loc, ms), distMedian(pr, loc, ms), distN(pr, loc, ms) );
                        mnStd_str = str_mn_std_n(distMean(pr, loc, ms), distStderr(pr, loc, ms), N_ifShow );
                        
                        if ~isempty(pval_str) && showSigTestsOnPlots
                            title_str = {mnStd_str, pval_str};
                        else
                            title_str = mnStd_str;
                        end
%                                      str_mn_std2(distMean(pr, loc, ms), distStderr(pr, loc, ms), round(distN(pr, loc, ms)) );
                        tun_flash_dphi = strcmp(cmpType_s, 'phase') && strcmp(gratingType_s, 'flashed') && any(strcmp(measures{ms}, {'dphi', 'dF1'}));
                        if histFirstTime || ~h_ax(pr, loc, ms)
                            figure(pr+gratingType*10);
                            h_ax(pr, loc, ms) = subplotGap(subM, subN, m_i, loc_i, M_spacing, N_spacing); hold on;
                            
                            binEdges = S{loc, ms}.binEdges;                                                        
                            binCenters = binEdge2cent(binEdges);
                            
                            
                            barW = iff(tun_flash_dphi, .65, 1);
%                             barW = 1;
                            h_bars(pr, loc, ms, 1:nMaxCats) = bar(binCenters, bin_N{pr, loc, ms}, barW, 'barlayout', 'stacked');
                            if ~isempty(bin_dN{pr, loc, ms})
                                hold on;
                                h_err(pr, loc, ms) = errorbar(binCenters, sum(bin_N{pr, loc, ms},2), bin_dN{pr, loc, ms}, '.');
                                hold off
                            end

                            xlim(h_ax(pr, loc, ms), [binEdges(1), binEdges(end)]);
                            if strcmp(cmpType_s, 'phase') && any(strcmp( measures{ms}, {'dphi', 'dF1'}) ) 
                                set(gca, 'xtick', iff(tun_flash_dphi && (length(binCenters) < 15), 0:45:180, binCenters) );
                            end                            
                            colormap(hsv(nMaxCats+1));
                        else                            
                            if tun_flash_dphi
                                binEdges = S{loc, ms}.binEdges;                                                        
                                binCenters = binEdge2cent(binEdges);
                                if (length(binCenters) < 15)
                                    xticks_dphi =  0:45:180;
                                    barWidth = 0.8;
                                else
                                    xticks_dphi =  0:30:180;
                                    barWidth = 1;
                                end
                                set(h_ax(pr, loc, ms), 'xtick', xticks_dphi );
                                set(h_bars(pr, loc, ms, 1:nMaxCats), 'barwidth', barWidth);
                            end
%                                 set(h_ax(pr, loc, ms), 'xtick', iff((length(binCenters) < 15), 0:45:180, binCenters) );
                                
                            adjustStackedBarHeights(h_bars(pr, loc, ms, 1:nMaxCats), bin_N{pr, loc, ms});
                            if ~isempty(bin_dN{pr, loc, ms})
                                set(h_err(pr, loc, ms), 'ydata', sum(bin_N{pr, loc, ms},2), 'udata', bin_dN{pr, loc, ms}, 'ldata', bin_dN{pr, loc, ms});
                            end

                            if adjustAxesLocation
                                subplotGap(subM, subN, m_i, loc_i, M_spacing, N_spacing, h_ax(pr, loc, ms));
                            end
                            
                            barCol = switchh(pairTypes{pr}, {'Wcc', 'Bcc'}, ['b', 'r', 'm']);
                                
                            set(h_bars(pr, loc, ms, 1:nMaxCats), 'facecolor', barCol)
                            
                        end
                        fontsize = 11;
                        h_title(pr, loc, ms) = updateTitle(h_title(pr, loc, ms), title_str, 'fontsize', fontsize+1, 'parent', h_ax(pr, loc, ms));

                        if (loc_i == 1) && ~skipYlabels
                            h_ylabel(pr, loc, ms) = axlabel('y', h_ylabel(pr, loc, ms), h_ax(pr, loc, ms), measureLabels{ms});  
                            set(h_ylabel(pr, loc, ms), 'rotation', 0, 'horizontalAlignment', 'right', 'visible', 'on');
                        end
                        if skipYlabels
                            set(h_ylabel(pr, loc, ms), 'visible', 'off');
                        end
                        
                        if m_i == length(measureInds) && (showLocLabelsIfOneLabel || (nLoc > 1))
                            h_xlabel(pr, loc, ms) = axlabel('x', h_xlabel(pr, loc, ms), h_ax(pr, loc, ms), locationLabels{loc}); 
                        end
                        
                        % fix y limits of all axes of this measure:
                        equalizeMeasureYLims = 0;
                        if equalizeMeasureYLims
                            ylim_fields = {'YLim', m_ylim2(ms,:)};
                        else
                            ylim_fields = {'YLimMode', 'auto'};
                        end
                        set(h_ax(pr, loc, ms), ylim_fields{:}, 'YTickMode', 'auto', 'YTickLabelMode', 'auto', 'tickdir', 'out', 'box', 'off');
                        yticklabels = num2str(get(h_ax(pr, loc, ms),'YTick').'); %% skip the "x10^n" auto scaling
                        set(h_ax(pr, loc, ms),'YTickLabel', yticklabels);
                            
                    end
                                            
%                             if (loc_i == 1)
%                                 if (h_ylabel(pr, loc, ms) == 0)
%                                     h_ylabel(pr, loc, ms) = ylabel(h_ax(pr, loc, ms), measureLabels{ms});
%                                 end
%                             end
%                                                         
%                             if (m_i == length(measureInds))
%                                 if (h_xlabel(pr, loc, ms) == 0),
%                                     h_xlabel(pr, loc, ms) = xlabel(h_ax(pr, loc, ms), locationLabels{loc});
%                                 end
%                             end


                        
                end
            
                
                
                % things done once per figure:
                % Set Main Figure Title                                
                if ~skipBigHeading
                    showNinTitle = false;
                    if showNinTitle
                        allNs = distN(pr, :, :); allNs = allNs(activeLocMs); maxN = max(allNs); 
                        Ntext = ['  (N = ' num2str(maxN) ')'];
                    else
                        Ntext = '';
                    end                                
                    titleStr = [titleCase( gratingType_s ) ' Gratings ' pairTypeLabels{pr}  Ntext];
                    if histFirstTime                    
                        h_figtitle(pr) = suptitle_2( titleStr );
                    else
                        set(h_figtitle(pr), 'string', titleStr);                    
                    end
                    set(h_figtitle(pr), 'fontsize', fontsize+5);
                end
                
                % Set / Update Legend                
                if ~hideLegend
                    testBars = false;                     
                    p_axUR = get(h_ax(pr, locationInds(end), measureInds(1) ), 'position');
                    dummy_axpos = iff(testBars, [.9, .5, .1, .2], [p_axUR(1)+p_axUR(3), p_axUR(2)+p_axUR(4), .01, .01]);
                    curLegend = colorCodings(cc_i).labels{pr};
                    if histFirstTime || ~ishandle(h_dummy_ax(pr))
                        figure(pr+gratingType*10);
                        h_dummy_ax(pr) = axes('position', dummy_axpos);                
                        h_dummy_bars(pr,:) = bar([1 2], ones(2, nMaxCats));
                        
                        idxNonempty = cellfun(@(x) ~isempty(x), curLegend);                      
                        h_legend(pr) = legend(h_dummy_bars(pr,idxNonempty), curLegend(idxNonempty), 'location', 'NorthWest', 'fontsize', legendFontSize);
                        set([h_dummy_ax(pr) h_dummy_bars(pr,:)], 'visible', 'off'); 
                    end               
                    if updateLegend
                        set(h_dummy_bars(pr,:), 'visible', 'on');
                        if fastLegendUpdates
                            set(h_legend(pr), 'String', curLegends);
                        else
                            idxNonempty = cellfun(@(x) ~isempty(x), curLegend);                      
                            h_legend(pr) = legend(h_dummy_bars(pr,idxNonempty), curLegend(idxNonempty), 'location', 'NorthWest', 'fontsize', legendFontSize);                    
                        end

    %                     set([h_dummy_ax(pr) h_dummy_bars(pr,:)], 'visible', 'off');
                        set(h_dummy_ax(pr), 'visible', 'off');
                        set(h_dummy_bars(pr,:), 'visible', 'off');
                    end
                    leg_pos = get(h_legend(pr), 'position');                
                    set(h_legend(pr), 'position', [p_axUR(1)+p_axUR(3)+.01, p_axUR(2)+p_axUR(4)-leg_pos(4), leg_pos(3:4)]);
                end
%                 if hideLegend
%                     set(h_legend(pr), 'visible', 'off');
%                 end
                    
                
                % align all axes
                for loc = locationInds                
                    alignSizes(h_ax(pr,loc,measuresActive), 'vertical')
                end
                for ms = measureInds
                    alignSizes(h_ax(pr,locationsActive,ms), 'horizontal')
                end
                
            end
            
            if ~hideLegend
                legendsVisible = arrayfun(isVisible, h_legend);
                legendsToShow = h_legend & ~legendsVisible & (pairTypesActive' & ~strcmp(colorCodingVariable, '[none]'));
                legendsToHide = h_legend & legendsVisible & (~pairTypesActive' | strcmp(colorCodingVariable, '[none]'));

                cellfun(@(hs) set(hs, 'visible', 'on'),  {h_legend(legendsToShow) });
                cellfun(@(hs) set(hs, 'visible', 'off'), {h_legend(legendsToHide) });                        
            end
                        
            histFirstTime = false;
            [prevSubM, prevSubN] = deal(subM, subN);
        end
        
        % summary figure;
        if showControlPanel && showSummaryFig            
        
            myNum2Str = @(x, d) iff(abs(x-round(x))<1e-10, num2str(round(x), '%d'), num2str(x, ['%.' num2str(d) 'f']) );         
            headingStr = [titleCase(gratingType_s) ' Gratings'];
            if nnz(filterActive) > 0
                acFilterNames = cellfun(@(s) strrep(s, '_', '\_'), filterFields(filterActive), 'un', 0);
                acFilterOps = cellfun(@(f) cmpFunc2Str(f), filterOps(filterActive), 'un', 0);
                acFilterTh = cellfun(@(s) myNum2Str(th.(s), 2), filterFields(filterActive), 'un', 0);
                blnk = repmat({' '}, 1, length(acFilterNames)); cma = [repmat({'  ||  '}, 1, length(acFilterNames)-1), ' '];
                C = [acFilterNames; blnk; acFilterOps; blnk; acFilterTh; cma];
                filterStr = [headingStr ' : ' C{:}];
            else
                filterStr = [headingStr];
            end
            
            y_top = .5;
%             x_lft = .5;
            boxDataW = .12;
            boxTotalW = .25;    

            if summaryFirstTime
                figure(sumFigId); clf; axis ij;
                sumFigAxes = gca;
                set(sumFigAxes, 'position', get(sumFigAxes, 'outerPosition'), 'xtick', [], 'ytick', []);
            else
                axes(sumFigAxes);
            end

            for p_i = 1:length(pairTypeInds)
                pr = pairTypeInds(p_i);
                
                for loc_i = 1:length(locationInds)            
                    loc = locationInds(loc_i);
                    for m_i = 1:length(measureInds)
                        ms = measureInds(m_i);

%                         mss = {distMean(pr, loc, ms), distStd(pr, loc, ms), distStderr(pr, loc, ms)};
                        mss = {distMean(pr, loc, ms), 1, distStderr(pr, loc, ms)};
                        
                        loc_offset = ((length(measureInds)+1)*(loc_i-1)*sumFigSpacing);
                        x_loc = loc_offset + sumFigSpacing*m_i;
%                         y_loc = find(strcmp(sumFigPairRows, pairTypes{pr})) + .5;
                        y_loc = p_i + y_top-.4;
                        
                        ax = [-1,  ((length(measureInds)+1)*length(locationInds)*sumFigSpacing), 0, nnz(pairTypesActive)+y_top];
                        locLab_x = 1+loc_offset + length(measureInds)*sumFigSpacing/2;
%                         msLab_x  = x(length(measureInds)*sumFigSpacing/2) + loc_i*length(measureInds)*sumFigSpacing;
                        
                        axis(sumFigAxes, ax);
                        sf_text = {str_mn_std_n(distMean(pr, loc, ms), distStderr(pr, loc, ms), distN(pr, loc, ms) ), ...
                                    makePvalStr(distPval{pr, loc, ms}, pvalLabels{pr,loc,ms}) };
                        h_sumFigBox{pr, loc, ms} = stdevBox(h_sumFigBox{pr, loc, ms}, x_loc, y_loc, mss, boxDataW, boxTotalW, pairColors{pr});
                        h_sfText(pr, loc, ms) = updateText(h_sfText(pr, loc, ms), x_loc, y_loc-.5, sf_text, ...
                            'fontsize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top' );
                        
                        if (loc_i == 1 && m_i == 1)
                            h_sfPairLabel(pr) = updateText(h_sfPairLabel(pr), 0, y_loc, pairTypes{pr}, 'HorizontalAlignment', 'right');
                        end
                        if p_i == 1 && m_i == 1
                            h_sfLocLabel(loc) = updateText(h_sfLocLabel(loc), locLab_x, .2, locationLabels{loc}, 'HorizontalAlignment', 'center');
                        end
                        if p_i == 1
                            h_sfMsLabel(loc, ms) = updateText(h_sfMsLabel(loc, ms), x_loc, .4, measureLabels{ms}, 'HorizontalAlignment', 'center');
                        end                        
                        if (loc_i == 1 && m_i == 1 && p_i == 1)
%                             yy = y_top -.1 +sum( cellfun(@(x) ~isempty(x),  strfind(pairTypes(pairTypesActive), 'cc')) );
%                             h_sfLine = updateLine(h_sfLine, 'xdata', [-10; 100], 'ydata', [yy;yy], 'color', 'k');                            
                            set(h_sfLine, 'visible', 'off')
                            h_sfFilterStr = updateText(h_sfFilterStr, mean(ax(1:2)), 0.02, filterStr, 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                        end
                            
                    end
                    
                end
            end
            summaryFirstTime = false;
                        
            cellfun(@(hs) set(hs, 'visible', 'on'), ...
                {h_sfPairLabel(pairTypesActive), h_sfLocLabel(locationsActive), h_sfMsLabel(activeLocMs), ...
                [h_sumFigBox{activePsLocMs}], h_sfText(activePsLocMs) });

            cellfun(@(hs) set(hs, 'visible', 'off'), ...
                {h_sfPairLabel(~pairTypesActive), h_sfLocLabel(~locationsActive), h_sfMsLabel(~activeLocMs), ...
                [h_sumFigBox{~activePsLocMs}],  h_sfText(~activePsLocMs)});

%             cellfun(@(hs) set(hs, 'visible', 'off'), {h_ax(~activePsLocMs), h_bars(~activePsLocMs), h_sumFigBox{~activePsMsLoc}
%                 h_ylabel(~activePsLocMs), h_xlabel(~activePsLocMs)},...

        end
        
        
%         distN
        3;
        imagePlotDphiDf1Corr = false;
        if imagePlotDphiDf1Corr
            
            idx_ms_x = find(strcmp(measures, 'dphi'));
            idx_ms_y = find(strcmp(measures, 'dF1'));
            
            for p_i = 1:length(pairTypeInds)
                pr = pairTypeInds(p_i);

                Xvals = setVals{pr, 1, idx_ms_x};
                Yvals = setVals{pr, 1, idx_ms_y};
                
                if isempty(Xvals) || isempty(Yvals)
                    distStderr(pr, 1, idx_ms_x) = nan;
                else                
                    distStderr(pr, 1, idx_ms_x) = corr(Xvals(:), Yvals(:));
                end
            end
        end
                
                
%       bin_N,    bin_dN,    distMean,    distStd,    distStderr,    distMedian,    distN,    distPval,    setVals,    setIdx
                        
        
        if showMsVsMs
%             idx_dphi = find(strcmp(measures, 'cc'));
%             idx_dF1 = find(strcmp(measures, 'rho'));
            idx_ms_x = find(strcmp(measures, ms_x));
            idx_ms_y = find(strcmp(measures, ms_y));
            idx_wcc = find(strcmp(pairTypes, 'Wcc'));
            if length([idx_wcc]) == 1 %ie none are empty
%                 X = setVals{idx_wcc, 1, idx_ms_x};
%                 Y = setVals{idx_wcc, 1, idx_ms_y};
                idx_ms_xy = intersect( setIdx{idx_wcc, 1, idx_ms_x}, setIdx{idx_wcc, 1, idx_ms_y});
                X = S{idx_ms_x}.val(idx_ms_xy);
                Y = S{idx_ms_y}.val(idx_ms_xy);

                corr_r = corr(X(:), Y(:), 'type', 'pearson');
%                 corr_rho = corr(X(:), Y(:), 'type', 'spearman');
                
                msVsMsFirstTime = true;
                if msVsMsFirstTime
                    
                    figure(msVsMsFig); clf; 
%                     h_ms12_ax = subplot(2,1,1);
                    h_ms12_L = plot(X, Y, 'bo', 'markersize', 2);
                    
                    xlabel(measureLabels{idx_ms_x});         
                    ylabel(measureLabels{idx_ms_y});
                    xlims = S{idx_ms_x}.binEdges([1, end]);  
                    ylims = S{idx_ms_y}.binEdges([1, end]);                    
                    xlim(xlims + [-1,1]*diff(xlims)/30);     
                    ylim(ylims + [-1,1]*diff(ylims)/30);
                    
                    xticks = iff(xlims(2) > 90, [0:45:180], linspace(-1, 1, 11));
                    yticks = iff(ylims(2) > 90, [0:45:180], linspace(-1, 1, 11));
                    set(gca, 'xtick', xticks, 'ytick', yticks);

                    title(sprintf('r = %.2f (N = %d)', corr_r, length(X)));
                    set(msVsMsFig, 'windowButtonDownFcn', {@selectPointsInFigure, @updateTuningCurvePair});

                    figure(msVsMsFig+1); clf; 
%                     h_ms_ex_ax = subplot(2,1,2);
                    h_ms12_L = plot(0,0, 'bo-', 0,0, 'bs--',  0,0, 'go-', 0,0, 'b:', 0,0, 'g:', 0,0, 'bd', 0,0, 'gd');
                    set(h_ms12_L(1), 'linewidth', 2, 'markerfacecolor', 'b')
                    set(h_ms12_L(3), 'linewidth', 2, 'markerfacecolor', 'g')
                    set(h_ms12_L(6:7), 'markerfacecolor', 'auto');
                    h_ms12_tit = title(' ');
                    xlim([0 360]);
                    set(gca, 'xtick', [0:90:360]);
                else
                    set(h_ms12_L, 'xdata', X, 'ydata', Y);    
                    set(get(h_ms12_L, 'parent'), 'ylim', [-1 1]);
                    
                end                
                
            end    
            
        end
    

        firstTime = false;
        toDisplayPrev = toDisplayCur;

        colorCodingVariable_prev = colorCodingVariable;
        [curBin_N, curBin_dN, curDistMean, curDistStd, curDistStderr, curDistMedian, curDistN, curDistPval, curSetVals, curSetIdx] = deal(...
         bin_N,    bin_dN,    distMean,    distStd,    distStderr,    distMedian,    distN,    distPval,    setVals,    setIdx);
        curPvalLabels = pvalLabels;
     
        setFilterActive_prev(activePsLocMs) = {filterActive};
        setTh_prev(activePsLocMs) = {th};
        nResamplesPrev = nResamples;        
    end


    function computeAndPrintSigTests(~,~)
        
        
        
        
        
        
        
    end    
        
       

    function updateTuningCurvePair(glob_id, grp_id, loc_id, pt_ms_x, pt_ms_y) %#ok<INUSL>
        
%         Gid_find = 1153;
%         cellIds_find = [2,5];
%         allGroupIds = cat(1, pairData.Gids);
%         allCellIds = cat(1, pairData.cellIds);
%         pair_idx = find(allGroupIds(:,1) == Gid_find & allGroupIds(:,2) == Gid_find & allCellIds(:,1) == cellIds_find(1) & allCellIds(:,2) == cellIds_find(2) );
        
        
        pair_idx = idx_ms_xy(loc_id);%curSetIdx{idx_wcc, 1, idx_ms_x}(loc_id);
        pData = pairData( pair_idx );        
        assert(pt_ms_x == S{idx_ms_x}.val(pair_idx));
        assert(pt_ms_y == S{idx_ms_y}.val(pair_idx));

        Gid = pData.Gids(1);
        cellIds = pData.cellIds;
        
        tc1 = pData.phaseTCs(:,1,1); tc1 = tc1/max(tc1);
        tc2 = pData.phaseTCs(:,1,2); tc2 = tc2/max(tc2);
        nPh = length(tc1);
        ph_ext = linspace(0,360, nPh+1);  ph = ph_ext(1:nPh);
        dph = round(360/nPh);
                
        dphi_bin = deltaPhi(ph, tc1, tc2, 'cross-correlation') / dph;
        if dphi_bin < 0
            dphi_bin = nPh+dphi_bin;
        end
        tc1_shft = shft(tc1, dphi_bin);
                
        [phi1_rad, f_cos1, t_f1] = getF1phase( ph, tc1, 360);
        [phi2_rad, f_cos2, t_f2] = getF1phase( ph, tc2, 360);
        phi1_deg = rad2deg(phi1_rad); phi1_y = f_cos1(indmin(abs(t_f1-phi1_deg)));
        phi2_deg = rad2deg(phi2_rad); phi2_y = f_cos2(indmin(abs(t_f2-phi2_deg)));
                
        f1_show_tf = isequal(sort([idx_ms_x, idx_ms_y]), [3 4]); %'off';
        f1_show = iff(f1_show_tf, 'on', 'off');
        tc1_shift_show = 'off';
        set(h_ms12_L(1), 'xdata', ph_ext, 'ydata', wrp(tc1));
        set(h_ms12_L(3), 'xdata', ph_ext, 'ydata', wrp(tc2));        
        set(h_ms12_L(2), 'xdata', ph_ext, 'ydata', wrp(tc1_shft), 'visible', tc1_shift_show);
        set(h_ms12_L(4), 'xdata', t_f1, 'ydata', f_cos1, 'visible', f1_show);
        set(h_ms12_L(5), 'xdata', t_f2, 'ydata', f_cos2, 'visible', f1_show);
        set(h_ms12_L(6), 'xdata', phi1_deg, 'ydata', phi1_y, 'visible', f1_show);
        set(h_ms12_L(7), 'xdata', phi2_deg, 'ydata', phi2_y, 'visible', f1_show);
        set(h_ms12_L(6), 'markerfacecolor', 'auto', 'marker', '*', 'color', 'r');%[.2 .2 1]);        
        set(h_ms12_L(7), 'markerfacecolor', 'auto', 'marker', '*', 'color', 'r');%[.2 1 .2]);        
        ax = get(h_ms12_tit, 'parent');
        legend off;
%         legend(ax, h_ms12_L([1,3]), {['Cell ' num2str(cellIds(1))], ['Cell ' num2str(cellIds(2))]}, 'location', 'best')
%         tc_cc = pearsonR(tc1, tc2);
%         tc_rho = spearmanRho(tc1, tc2);
%         tc_dphi = deltaPhi(ph, tc1, tc2, 'cross-correlation');
%         tc_dF1 = deltaPhi(ph, tc1, tc2, 'angle');
        [cc, rho, dp, df1] = deal(S{1}.val(pair_idx), S{2}.val(pair_idx), S{3}.val(pair_idx), S{4}.val(pair_idx) );
        s1 = sprintf('Site # %d. Cells %d, %d.', Gid, cellIds(1), cellIds(2));        
        ms_strs = {sprintf('cc = %.2f', cc), sprintf('\\rho = %.2f', rho), sprintf('\\Delta\\phi = %.0f', dp), sprintf('\\DeltaF1 = %.0f', df1)};
        s2 = [ms_strs{idx_ms_x} ', ' ms_strs{idx_ms_y}];
%         s2 = sprintf('r = %.2f, \\rho = %.2f, \\Delta\\phi = %.0f, \\DeltaF1 = %.0f', cc, rho, dp, df1);
        set(h_ms12_tit, 'string', {s1, s2}, 'fontsize', 14);%[.2 1 .2]);        
        
        
    end

    function [distMean, distStd, distStderr, distN, distPval, bin_N, bin_dN, setVals, setIdx] = updateFilters(varargin)
        
        [toDisplayCur, ms_x, ms_y, showPairTypes, nResamples, showLocations, showMeasures, colorCoding] = varargin{1:8};
        filterVarArgs = varargin(9:end);
        
        [showHistograms, showSummaryFig, showMsVsMs] = dealV( toDisplayCur );
        
        pairTypesActive = showPairTypes;        
        locationsActive = showLocations;
        measuresActive  = showMeasures;
        colorCodingVariable = colorCoding;
        
        %             pairFilters1,    minRepPval_cc_p, nPhases,  ...
%             pairFilters2,   SCtypeStr, minF1oDC_pref_lo, minF1oDC_pref_hi, ...
%             pairFilters3,   animalCmpStr, penetrCmpStr, locCmpStr, ...
%             locFilters1,    minF1oDC_cmp_lo, minF1oDC_cmp_hi, minSumPhs, ...
%             locFilters2,    minFracR_lo, minFracR_hi, ...
%             setFilters,     pval] = varargin{:};
        
        filterActive = cat(1, filterVarArgs{filterHeadersIdx})';
        
%             [limByRepPval_av, limByNPhases] = dealV(pairFilters1);
%             [limBySCtype, limByF1oDC_pref_lo, limByF1oDC_pref_hi] = dealV(pairFilters2);
%             [limByAnimalCmp, limByPenetrCmp, limBySameLoc] = dealV(pairFilters3);
%             [limByF1oDCs_cmp_lo, limByF1oDCs_cmp_hi, limBySumPhs] = dealV(locFilters1);
%             [limByMinFracR_lo, limByMinFracR_hi] = dealV(locFilters2);
%             [limByPval] = dealV(setFilters);
              
                
        
%         filterActive = [limBySCtype, limByF1oDC_pref_lo, limByF1oDC_pref_hi, limByRepPval_av, limByNPhases, limByAnimalCmp, limByPenetrCmp, limBySameLoc ...
%                         limByF1oDCs_cmp_lo, limByF1oDCs_cmp_hi, limBySumPhs, limByMinFracR_lo, limByMinFracR_hi, ...
%                         limByPval ];
        
%         SCtype = length( strfind(SCtypeStr, 'S') );  
%         
%         animalCmp = switchh(animalCmpStr, {'same', 'diff'}, [1, 0]);
%         penetrCmp = switchh(penetrCmpStr, {'same', 'diff'}, [1, 0]);
%         locCmp    = switchh(locCmpStr,    {'same', 'diff'}, [1, 0]);
        
%         animalCmpPairId = find(strcmp(pairFilterVars, 'animalCmp'));
%         switch animalCmpStr
%             case 'diff', animalCmp = 0.5+rand/1000; pairFilterVarOp{animalCmpPairId} = @le;
%             case 'same', animalCmp = 0.5+rand/1000; pairFilterVarOp{animalCmpPairId} = @ge;
% %             case 'auto', animalCmp = 1.5+rand/1000; pairFilterVarOp{animalCmpPairId} = @gt;
%         end        

        % modify values for string arguments (that must be converted to numeric values 
        strArgIdxs = find(cellfun(@ischar, filterVarArgs));
        for str_i = strArgIdxs;
            strVal = filterVarArgs{str_i};
            switch filterArgNames{str_i}
                case 'n_phases', numVal = str2double(strVal);   %'4' -> 4, '8' -> 8;
                case {'SCtype_pref', 'SCtype_cmp'},   numVal = length( strfind(strVal, 'S') );  % 'C-C' -> 0,  'S-C'-> 1, 'S-S' -> 2.
                case {'animalCmp', 'penetrCmp', 'locCmp'}, numVal = switchh(strVal, {'same', 'diff'}, [1, 0]);            
            end
            filterVarArgs{str_i} = numVal;
        end
        
%         % modify values for n_phases type filter:
%         nPhases_idx = find(strcmp('n_phases', filterArgNames));
%         nPhaseStr = filterVarArgs{SCtype_idx};
%         SCtype = length( strfind(SCtypeStr, 'S') );
%         filterVarArgs{SCtype_idx} = SCtype;
% 
%         % modify values for S/C type filter:
%         SCtype_idx = find(strcmp('SCtype', filterArgNames));
%         SCtypeStr = filterVarArgs{SCtype_idx};
%         SCtype = ;
%         filterVarArgs{SCtype_idx} = SCtype;
%         
%         % modify values for animalCmp, penetrationCmp, locationCmp filters:
%         locsCmp = strfind(filterArgNames, 'Cmp');
%         for cmp_idx = find(cellfun(@(x) ~isempty(x), locsCmp));
%             cmpStr = filterVarArgs{cmp_idx};
%             cmpVal = switchh(cmpStr, {'same', 'diff'}, [1, 0]);
%             filterVarArgs{cmp_idx} = cmpVal;
%         end

        thStructArgs =  [filterArgNames(filterValsIdx); filterVarArgs(filterValsIdx)];
        th = struct(thStructArgs{:});

%         th = struct('SCtype', SCtype, 'minF1oDC_pref_lo', minF1oDC_pref_lo,  'minF1oDC_pref_hi', minF1oDC_pref_hi, 'min_rep_cc_p', minRepPval_cc_p, ...
%             'n_phases', str2double(nPhases), 'animalCmp', animalCmp, 'penetrCmp', penetrCmp, 'locCmp', locCmp, ...
%                     'minF1oDC_cmp_lo', minF1oDC_cmp_lo, 'minF1oDC_cmp_hi', minF1oDC_cmp_hi, 'minSumPhs', minSumPhs,  'minFracOfMaxes_lo', minFracR_lo, 'minFracOfMaxes_hi', minFracR_hi,...
%                      'pval', pval );
        th = structfun(@double, th, 'un', 0);
        assert( isequal(fieldnames(th), filterFields') ); % make sure are in correct order        
        
        [distMean, distStd, distStderr, distN, distPval, bin_N, bin_dN, setVals, setIdx] = plotOrUpdateAllFigures;
%         fprintf('Finished replotting: '); toc;
        
    end    
    
    
    tf = [true, false];    
%     pvals = [0 : .01 : 1];    
%     pvals = [0:.01:1];    
%     strvals = [-1:.01:1];   
    
    
    showPairOptions = repmat({tf}, 1, nP);   showPair0 = true(1,nP);
    showLocOptions  = repmat({tf}, 1, nLoc); showLoc0  = true(1,nLoc);
    showMsOptions   = repmat({tf}, 1, nM);   showMs0   = true(1,nM);
    
%     rng02 = linspace(0:.05:2);    
        
    
    pairTypeDepVars = cell(1,nP);
    pairTypeDepVars(:) = {{{},{}}};
    pairTypeDepVars( (~cellfun(@isempty, strfind(pairTypes, '_r'))) ) = {{{}, {'nResamples'}}};
        
    func_handle = @updateFilters;
    tfDepVars = @(varname) { {}, {varname} };
    

    % Generate Filter-related Arguments for "manipulate"
    
    [filterTypes, idx_orig] = unique({allFilters.type});
    filterTypes = filterTypes(ord(idx_orig)); % put back in original order.
    nFilterArgs = length(allFilters)+length(filterTypes);
    f_arg_i = 1;
    filterArgs = cell(1,nFilterArgs);    
    for ft_i = 1:length(filterTypes)
        ft_idx = find( strcmp(filterTypes{ft_i}, {allFilters.type}) );
        ft_grps = [allFilters(ft_idx).group];
        uFt_grps = unique(ft_grps);
        for grp_id = uFt_grps            
            grp_filter_idxs = find(ft_grps == grp_id);
            numFiltersInGrp = length(grp_filter_idxs);
            mainLabel = [filterTypes{ft_i}, 'Filter' num2str(grp_id)];            
            mainTF = repmat(tf, numFiltersInGrp, 1);
            mainInitF = false(numFiltersInGrp, 1);
            mainDepVars = cellfun(tfDepVars, {allFilters(ft_idx(grp_filter_idxs)).filterName}, 'un', 0);
%             if numFiltersInGrp == 1
%                 mainDepVars = mainDepVars{1};                
%             end
            mainLabels  = {allFilters(ft_idx(grp_filter_idxs)).label};
            filterArgs{f_arg_i} = {mainLabel, mainTF, mainInitF, mainDepVars, mainLabels};
            f_arg_i = f_arg_i+1;
            for f_i = grp_filter_idxs
                fltr = allFilters(ft_idx(f_i));
                if iscellstr(fltr.range)
                    filterArgs{f_arg_i} = {fltr.filterName, fltr.range, fltr.range{1}};
                else
                    filterArgs{f_arg_i} = {fltr.filterName, fltr.range, fltr.init};
                end
                f_arg_i = f_arg_i+1;
            end
        end
    end
    filterArgNames = cellfun(@(C) C{1}, filterArgs, 'un', 0);                            
    filterHeadersIdx = cellfun(@(s) ~isempty(strfind(s, 'Filter')), filterArgNames);        
    filterValsIdx = ~filterHeadersIdx;

%     figure(msVsMsFig);
%     h_mm = plot(0,0, 'bo');
    toDisplayOptions = repmat({[false,true]}, 1,3);
    toDisplayDepVars = { {{}, {}}, {{}, {}}, {{}, {'ms_x', 'ms_y'}} }; 
    
    args = { { 'toDisplay', toDisplayOptions, toDisplayCur, toDisplayDepVars, toDisplayNames }, ...
                 { 'ms_x', measures, measures{2} }, ...
                 { 'ms_y', measures, measures{2} }, ...
             { 'pairTypes', showPairOptions, showPair0, [], pairTypes, pairTypeDepVars}, ...
             { 'nResamples', [1:nMaxResamples], 1}, ...
             { 'locations', showLocOptions, showLoc0, [], locations }, ...
             { 'measures',  showMsOptions, showMs0, [], measures }, ...             
             {'colorCoding', {colorCodings.name} }, ...             
             filterArgs{:}, ...
            }; 
   
    if showControlPanel
        figure(panel_fig_id); clf;    
        set(panel_fig_id, 'Name', mainTitle);
%         manipulate(panel_fig_id, func_handle, args{:}, variableGroups, mainTitle);
        manipulate(func_handle, args, 'Title', mainTitle, 'FigId', panel_fig_id);
        
        h = uicontrol(panel_fig_id, 'style', 'pushbutton', 'units', 'pixel', 'position', [0 0 60 20], ...
            'string', 'Sig Tests', 'callback', @computeAndPrintSigTests);                
        
    end
    
        
%     function [distMean, distStd, distStderr, distN, distPval, bin_N] = updateFilters(inputs)

    function plotPairs(pairIdxs)
        allGids = [allCells.Gid];
        allCellIds = [allCells.cellId];
        wrp = @(x) x([1:end, 1]);
        ext = @(x) [x, x(end)+diff(x([1, 2]))];

        toPlot = 'psth';
        
        figure(13); clf;
        for i = 1:6
            pd = pairData(pairIdxs(i));
            Gid = pd.Gids(1);
            cellIds = pd.cellIds;
            cell1 = allCells( find(Gid == allGids & allCellIds == cellIds(1), 1) );
            cell2 = allCells( find(Gid == allGids & allCellIds == cellIds(2), 1) );
            subplot(3,2,i);

            switch toPlot
                case 'phase'
                    
                    ph_ext = ext(cell1.ph);                    
                    [ori_i, sp_i] = dealV(pd.oriSp_cmp);

                    r1 = squeeze( cell1.R(ori_i, sp_i, :) );
                    r2 = squeeze( cell2.R(ori_i, sp_i, :) );
                    plot( ph_ext, wrp(r1), 'bo-', ph_ext, wrp(r2), 'go-');
                    set(gca, 'xtick', ph_ext(1:2:end));
                    
                case 'psth';                   
                    [uori, usp, uph] = dbGetUniqueOriSpPh('Gid', Gid);
                    [nOri, nSp, nPh] = deal(length(uori), length(usp), length(uph));

                    binC = cell1.PSTH.bins;
                    [ori_i, sp_i] = dealV(pd.oriSp_cmp);
                    tst = zeros(nOri, nSp, nPh);
                    tst(ori_i, sp_i, :) = 1;
                    ori_sp_idx = find(tst(:));
                    
%                     vals1 = cell1.PSTH.vals;
%                     vals2 = cell2.PSTH.vals;
                    [bins, allVals1] = dbGetCellSpkStimHists(Gid, cellIds(1));
                    [bins, allVals2] = dbGetCellSpkStimHists(Gid, cellIds(2));                        
%                     new_order1 = ord(mean(allVals1, 1), 'descend');  
%                     new_order2 = ord(mean(allVals2, 1), 'descend');  
                    
                    vals1 = mean( allVals1(:,ori_sp_idx), 2);                    
                    vals2 = mean( allVals2(:,ori_sp_idx), 2);
                    
                    plot2PSTHs(bins, vals1, vals2);
                    binE = binCent2edge(binC);
%                     extt = @(x) [x, x(end)];
% 
%                     wh_pix = getObjDims(gca, 'pixels');
%                     x_range = binE(end)-binE(1);
%                     d = x_range / wh_pix(1);                    
%                     stairs(binE, y1, 'linewidth', 2); hold on; stairs(binE+d, y2, 'r', 'linewidth', 2);
                    xlim([binE(1), binE(end)]);
                    
            end
            id_str = sprintf('(%d) Gid = %d. Cells %d, %d', i, Gid, cellIds(1), cellIds(2));
            f1_dc_str = sprintf('F1/DCs: @max: %.1f, %.1f. @cmp: %.1f, %.1f', cell1.F1oDC_maxR_av, cell2.F1oDC_maxR_av, pd.loc_F1oDCs_cmp(1), pd.loc_F1oDCs_cmp(2));
            [cc, rho, dp, df1] = deal(S{1}.val(pairIdxs(i)), S{2}.val(pairIdxs(i)), S{3}.val(pairIdxs(i)), S{4}.val(pairIdxs(i)));
            val_str = sprintf('r = %.2f, \\rho = %.2f, \\Delta \\phi = %.0f, \\DeltaF1 = %.0f', cc, rho, dp, df1);
            title({id_str, f1_dc_str, val_str});        
            
        end
    end


    if makeParameterPlots || makeBarGraphOfResults
        if length(pairTypes) > 1
%             error('Warning: will take much longer if have Bcc active');
        end
        
        addSummaryDataOnFirstPlot = false;
        addPairTypeLabelToEachPlot = true;
        addAxesLabelsOnEachPlot = false;
        showSigTextOnAxes = false;
        addMeasureLabels = true;
        saveDataToFile = true;
%         saveResults
        
        pairTypes_h = pairTypes;        
%         pairTypes1 = {'Wcc', 'Bcc'};
        locations_h = {'maxMinFracR'};

        measures_h = {'cc', 'rho', 'dphi', 'dF1'};
%         measures_h = {'cc', 'dphi'};
%         measures_h = {'cc'};
%         measures_h = {'dphi', 'dF1'};

        pvalNamesForMs = {'p_{U}', 'p_{U}', 'p_{X}', 'p_{KS'};        
%         pvalNamesForMs = {'p_{U}','p_{X}', };        
                
%         dataToImage = 'data';
        dataToImage = 'pval';
        

        pairTypesActive = arrayfun(@(pt) any(strcmp(pt, pairTypes_h)), pairTypes);
        locationsActive = arrayfun(@(lc) any(strcmp(lc, locations_h)), locations);
        measuresActive = arrayfun(@(ms) any(strcmp(ms, measures_h)), measures);        
        
        nP_h = length(pairTypes_h);
        nLoc_h = length(locations_h);
        nM_h = length(measures_h);
        
        pairTypeInds_h = find(pairTypesActive);
        locationInds_h = find(locationsActive);
        measureInds_h = find(measuresActive);        
%                 measureInds_h = cellfun( @(ms) find(strcmp(measures, ms)),  measures_h);           
        measureLabels_h = measureLabels(measureInds_h);
        pvalLabels_h = pvalLabels(pairTypeInds_h, locationInds_h, measureInds_h);
        
        
        assert(length(pvalNamesForMs) == length(measures_h));
%         ind_0  = arrayfun(@(ms0) find( strcmp(ms0, measures_h), 1), {'cc', 'rho'});
%         ind_90 = arrayfun(@(ms0) find( strcmp(ms0, measures_h), 1), {'dphi', 'dF1'});
        
        % labels for plots
%         FD = [upper(gratingType_s(1)) ' : '];
%         pTypeNames = {'within_c_c', 'within_c_c_R', 'between_c_c', 'within_c_mu', 'within_c_mu_R', 'between_c_mu'};
    end
    
    if makeParameterPlots 
   
        if gratingType == 1
            axlims(1) = struct('name', 'mnF1oDCs_cmp',   'thName', 'minF1oDC_cmp_lo', 'vals', [0:.1:2], 'labl', 'min F1/DC (@cmp)');
            axlims(2) = struct('name', 'mxF1oDCs_cmp',   'thName', 'maxF1oDC_cmp_hi', 'vals', [0:.1:2], 'labl', 'max F1/DC (@cmp)');
        elseif gratingType == 2
            axlims(1) = struct('name', 'mnF1oDCs_cmp',   'thName', 'minF1oDC_cmp_lo', 'vals', [0:.1:2], 'labl', 'min F1/DC (@cmp)');
            axlims(2) = struct('name', 'mxF1oDCs_cmp',   'thName', 'maxF1oDC_cmp_hi', 'vals', [0:.1:2], 'labl', 'max F1/DC (@cmp)');
        end                
        axlims(end+1) = struct('name', 'F1oDCs_cmp1',  'thName', 'minF1oDC_cmp_lo', 'vals', [0,.5,1],'labl', 'min F1/DC (@cmp)');
%         axlims(end+1) = struct('name', 'F1oDCs_cmp', 'thName', 'minF1oDC_cmp', 'vals', [0:.2:1.2], 'labl', 'min F1/DC (@cmp)');
        f1odc_w = .3;
        axlims(end+1) = struct('name', 'F1oDCs_cmp_Wmn',  'thName', {{'minF1oDC_cmp_lo', 'minF1oDC_cmp_hi'}}, 'vals', {{[0:.05:2-f1odc_w], [f1odc_w:.05:2]}},'labl', ['F1/DC (cmp) window (' num2str(f1odc_w, '%.2f') ')']);
        axlims(end+1) = struct('name', 'F1oDCs_cmp_Wmx',  'thName', {{'maxF1oDC_cmp_lo', 'maxF1oDC_cmp_hi'}}, 'vals', {{[0:.05:2-f1odc_w], [f1odc_w:.05:2]}},'labl', ['F1/DC (cmp) window (' num2str(f1odc_w, '%.2f') ')']);
        f1odc_w = .5;
        axlims(end+1) = struct('name', 'F1oDCs_cmp_Wmnmx',  'thName', {{'minF1oDC_cmp_lo', 'maxF1oDC_cmp_hi'}}, 'vals', {{[0:.05:2-f1odc_w], [f1odc_w:.05:2]}},'labl', ['F1/DC (cmp) window (' num2str(f1odc_w, '%.2f') ')']);
        f1odc_max_hi = iff(gratingType == 1, 1.7, 1.5);
        axlims(end+1) = struct('name', 'F1oDCs_cmp_ranges1',  'thName', {{'minF1oDC_cmp_lo', 'maxF1oDC_cmp_hi'}}, 'vals', {{[0, 0, 0.75], [2, 0.7, f1odc_max_hi]}},'labl', ['F1/DC']);
        
        axlims(end+1) = struct('name', 'F1oDCs_pref', 'thName', 'minF1oDC_pref', 'vals', [0:.25:1.2], 'labl', 'min F1/DC (@pref)');
        axlims(end+1) = struct('name', 'repPvals',    'thName', 'minRepPval_cc_p', 'vals', [0:.2:3],    'labl', '-log(rep pval)');
%         axlims(5) = struct('name', 'fracRs',      'thName', 'minFracR',      'vals', [0:.2:.9],  'labl', 'min FracR');
        maxMinFracR = iff(gratingType == 1, 0.95, 1);
        axlims(end+1) = struct('name', 'fracRs',       'thName', 'minFracR_lo', 'vals', [0, .5:.05:maxMinFracR], 'labl', 'min FracR (lower bound)');        
        axlims(end+1) = struct('name', 'fracRs1',      'thName', 'minFracR_lo',      'vals', [0,.65, .9],  'labl', 'min FracR');
        frac_w = .2;
        axlims(end+1) = struct('name', 'fracRs_window',  'thName', {{'minFracR_lo', 'minFracR_hi'}}, 'vals', {{[0:.05:1-frac_w],[frac_w:.05:1]}}, 'labl', ['FracR window (' num2str(frac_w, '%.2f') ')']);        

        maxPhTcCC = iff(gratingType == 1, 0.8, 0.5);        
        axlims(end+1) = struct('name', 'phaseTC_cc',    'thName', 'minPhTC_cc', 'vals', [0:.05:maxPhTcCC], 'labl', 'minPhaseTC cc');        
        maxPhTcDot = iff(gratingType == 1, 0.95, 0.7);
        axlims(end+1) = struct('name', 'phaseTC_dot',    'thName', 'minPhTC_dot', 'vals', [0:.05:maxPhTcDot], 'labl', 'minPhaseTC dot');        

        %         if gratingType == 1
%         elseif gratingType == 2
%             axlims(6) = struct('name', 'fracRs1',   'thName', 'minFracR',      'vals', [0,.65, .85],  'labl', 'min FracR');
%         end
        axlims(end+1) = struct('name', 'SCtype_cmp',    'thName', 'SCtype_cmp',     'vals', {{'CC', 'SC', 'SS'}},'labl', 'S/C type');
        getInd = @(nm) find( arrayfun(@(s) strcmp(s.name, nm), axlims) );
        
        % choose which variables to display on X (and Y) axis.
%         X = axlims(getInd('F1oDCs_cmp'));

%         X = axlims(getInd('fracRs'));
%         X = axlims(getInd('F1oDCs_cmp'));
%         X = axlims(getInd('fracRs_window'));
%         X = axlims(getInd('fracRs_window'));

%         Y = axlims(getInd('F1oDCs_cmp_Wmnmx'));
%         Y = axlims(getInd('mxF1oDCs_cmp'));
        

        X = axlims(getInd('mnF1oDCs_cmp'));
        Y = axlims(getInd('mxF1oDCs_cmp'));

        
%         Y = axlims(getInd('F1oDCs_cmp_ranges1'));
%         Y = axlims(getInd('fracRs'));
%         Y = axlims(getInd('F1oDCs_cmp'));
%         Y = axlims(getInd('fracRs1'));


        Z = axlims(getInd('fracRs'));
%         Z = axlims(getInd('fracRs_window'));
%         Z = [];
        
        if isempty(X)
            error('unknown X');
        end
        if isempty(Y)
            error('unknown Y');
        end

        xs = X.vals; 
        multipleXth = iscell(X.thName);
        multipleYth = iscell(Y.thName);
        multipleZth = ~isempty(Z) && iscell(Z.thName);
        
        if multipleXth
            xs_cent = mean( [xs{1}; xs{2}], 1);
        else
            xs_cent = xs;
        end
        nx = length(xs_cent);
    
        if ~isempty(Y)
            ys = Y.vals;            
            if multipleYth
                ys_cent = mean( [ys{1}; ys{2}], 1);
            else
                ys_cent = ys;
            end            
            ny = length(ys_cent);
        else
            ny = 1;
        end

        if ~isempty(Z)
            zs = Z.vals;       
            if multipleZth
                nz = length(zs{1});
            else
                nz = length(zs);
            end
        else
            Z = struct([]);
            nz = 1;
        end
        
        
        multipleYvals = ny > 1;
        multipleZvals = nz > 1;
        
        doErrorBarPlots = ny <= 5 && nz == 1;
        doImagePlots = ny > 5 && nz == 1;
                
        
%         pairTypeIdxsHere = cellfun(@(s) find(strcmp(s, pairTypes)), pairTypesHere);
        
        % initialize: receiving variables.
        [data, data_std, numPairs, pvals] = deal( cell(nP_h, nLoc_h, nM_h) );
        data(:) = {zeros(nx, ny, nz)};
        data_std(:) = {zeros(nx, ny, nz)};
        numPairs(:) = {zeros(nx, ny, nz)};
        pvals(:) = {zeros(nx, ny, nz, nMaxSigTests)};     
        setVals = cell(nP_h, nx, ny, nz);
        setIdx  = cell(nP_h, nx, ny, nz);
%         [Wcc_maxR1xR2_cc,    Wcc_maxR1xR2_rho,    Wcm_maxR1xR2_cc,    Wcm_maxR1xR2_rho   ] = deal(;
%         [Wcc_maxMinFracR_cc, Wcc_maxMinFracR_rho, Wcm_maxMinFracR_cc, Wcm_maxMinFracR_rho] = deal(zeros(ny, nx));
%         [Wcc_NumPairs, Wcc_ccPval, Wcm_NumPairs, Wcm_ccPval] = deal(zeros(ny, nx));
        
%         t = true; f = false;        
        [toDisplay,nResamples, colorCoding] = deal([0 0 0], 1, '[none]');         

        basicInputs = {toDisplay, measures{1}, measures{2}, pairTypesActive, nResamples, locationsActive, measuresActive, colorCoding};
        
%         getDepVars = @(depList) cellfun(@(C) C{2}, depList);
%         depVars = cellfun(@(C) getDepVars(C{4}), filterArgs(filterHeadersIdx), 'un', 0)
        
%         filterArgs
%         filterArgNames 
%         filterHeadersIdx 
%         filterValsIdx ;
        filterInitVals = cellfun(@(C) C{3}, filterArgs, 'un', 0);
        
        filterInputs = cell2struct(filterInitVals, filterArgNames, 2);
        
        % initialize: set values        
%         filterInputs = enableFilter(filterInputs, 'min_rep_cc_p');
%         filterInputs.minRepPval_cc_p = 2;

%         filterInputs = enableFilter(filterInputs, 'minFracR_lo');
%         filterInputs.minFracR_lo = iff(gratingType == 1, .65, .7);
        
        
        
        % after initial compilation  set X & Y filters on.                        
        filterInputs = enableFilter(filterInputs, X.thName);
        if multipleYvals
            filterInputs = enableFilter(filterInputs, Y.thName);
        end
        if multipleZvals
            filterInputs = enableFilter(filterInputs, Z.thName);
        end
        
                
        
        % 1. CALCULATE
        progressBar('init-', nx*ny*nz); 
        for x_i = 1:nx
            filterInputs = setFilterValue(filterInputs, X.thName, xs, x_i);                    
            for y_i = 1:ny
                if multipleYvals,
                    filterInputs = setFilterValue(filterInputs, Y.thName, ys, y_i);
                end                
                for z_i = 1:nz
                    progressBar;  
                    
                    if multipleZvals,
                        filterInputs = setFilterValue(filterInputs, Z.thName, zs, z_i);
                    end
                    
                    filterInputs_C = struct2cell(filterInputs);
                    [distMean, distStd, distStderr, distN, distP, tmpBinN, tmpBinDN, sVals, sIdx] = updateFilters(basicInputs{:}, filterInputs_C{:});
                    
                    for p_j = 1:length(pairTypeInds_h)
                        setVals{p_j, x_i, y_i, z_i} = sVals(p_j,:,:);
                        setIdx{p_j, x_i, y_i, z_i} = sIdx(p_j,:,:);                    

                        for loc_j = 1:length(locationInds_h)
                            for m_j = 1:length(measureInds_h)
                                data{p_j, loc_j, m_j}(x_i, y_i, z_i)     = distMean(pairTypeInds_h(p_j), locationInds_h(loc_j), measureInds_h(m_j));
                                data_std{p_j, loc_j, m_j}(x_i, y_i, z_i) = distStderr(pairTypeInds_h(p_j), locationInds_h(loc_j), measureInds_h(m_j));
                                numPairs{p_j, loc_j, m_j}(x_i, y_i, z_i) = distN(pairTypeInds_h(p_j), locationInds_h(loc_j), measureInds_h(m_j));
                                pvals{p_j, loc_j, m_j}(x_i, y_i, z_i, :) = distP{pairTypeInds_h(p_j), locationInds_h(loc_j), measureInds_h(m_j)};
                            end

                        end
                    end                

                end
            end
        end
        progressBar('done');

        
        % 2. PLOT    
        if doImagePlots  % 2D plot - X on one axis, Y on second axis.
                                
%             if showPvalPlot
%                 n_ac = 2;
%                 nSubInds = {[3,1], [3,1]};
%                 pSubInds = {[3,2], [3,2]};            
%             else
%                 n_ac = 8;
%                 nSubInds = {[3,3], [3,6]};
%             end            
            jet1 = [1 1 1; jet(100)];
            minNPairs = 10;

            plotsOnBottom = {'NPairs'};
%             plotsOnBottom = {};
            anyPlotsOnBottom = ~isempty(plotsOnBottom);
            plot2ColsIfOneLoc = false;
            
            addAxesLabelsOnEachPlot = isempty(plotsOnBottom);            
            
            if plot2ColsIfOneLoc && (nLoc_h == 1) && ~odd(nM_h)
                M = floor(sqrt(nM_h));
                N = nM_h/M;
                M = M+anyPlotsOnBottom;
            else
                M = nM_h+anyPlotsOnBottom;
                N = nLoc_h;                
            end
            switch length(plotsOnBottom)
                case 1, nBottomPlots_idx = {subplotInds(M,4, [M,2], [M,3])};
%                 case 1, nBottomPlots_idx = {subplotInds(M,4, [M,1], [M,4])};
                case 2, nBottomPlots_idx = {subplotInds(M,4, [M,1], [M,2]),  subplotInds(M,4, [M,3], [M,4])};
%                 case 3, nBottomPlots_idx = {subplotInds(M,5, [M,1], [M,1]),  subplotInds(M,4, [M,3], [M,4])};
            end
            
                        
            h_ax = zeros(nP_h, nLoc_h, nM_h);
            for p_j = 1:nP_h;
                plot_j = 1;    
%                 allDataThisPair = [data{p_j,:,:}];
%                 clims = [0, max(allDataThisPair(:)) ];            
                
                p_s = iff(addPairTypeLabelToEachPlot, '', [pairTypes_h{p_j} ' : ']);
%                 pvals_plot = pvals{p_j,1,1}; pvals_plot(pvals_plot > .05) = nan;
                
                idx_ignore = numPairs{p_j,1,1} < minNPairs;
                
                figure(100+p_j+gratingType*10); clf; set(gcf, 'Name', [gratingType_s ' : ' pairTypes_h{p_j}]);
                for loc_j = 1:nLoc_h
                    for m_j = 1:nM_h                        
                        h_ax(plot_j) = subplot(M,N,plot_j);
%                         
                        
                        switch dataToImage
                            case 'pval', 
                                
                                pvalLabels_here = pvalLabels_h{p_j, loc_j, m_j};
                                pvalId = find(strcmp(pvalLabels_here, pvalNamesForMs{m_j}), 1);
                                if isempty(pvalId)
                                    pvals_here = nan(length(xs), length(ys));
                                else
                                    pvals_here = pvals{p_j, loc_j, m_j}(:,:,1,pvalId);
                                end
                                plotdata = -log10( pvals_here );
                                
                            case 'data', plotdata =data{p_j, loc_j, m_j}(:,:,1);
                        end
%                         plotdata = data_std{p_j, loc_j, m_j};
                     
                        if ~isempty(minNPairs)
                            plotdata(idx_ignore) = nan;                            
                        end
                                                
                        h_im(plot_j) = imagesc(xs_cent, ys_cent, plotdata'); % image the data with the pixels that we want blanked == nan.
%                         plotdata = adjustImageNans(plotdata);       
                        
                        axis xy;
                        colormap(jet1);
%                         title([p_s ' ' locations1{loc_j} ' : ' measures1{m_j} ]);                        
                        if (plot_j == 1) && (addSummaryDataOnFirstPlot)
                            title( [titleCase(gratingType_s) ' Gratings. frac > ' num2str(filterInputs.minFracR_lo, '%.2f')]);
                        end                        
%                         if (nLoc > 1) || showLocLabelsIfOneLabel                                                    
%                         end
                            

                        if addMeasureLabels
                            if addAxesLabelsOnEachPlot
                                title(measureLabels_h{m_j}) ;
                            else
                                ylabel(measureLabels_h{m_j}) ;
                            end
                        end
                        if addAxesLabelsOnEachPlot
                            ylabel(Y.labl); 
                            xlabel(X.labl); 
                        end
                        set(h_ax(plot_j), 'active', 'outer');
                            
                        h_col(plot_j) = colorbar;
%                         
                        plot_j = plot_j+1;
                    end
                end

                if ~plot2ColsIfOneLoc  % only do this when all on one column - for poster.
                    for loc_j = 1:nLoc_h
                        for m_j = 1:nM_h                        
                            p_ax = get(h_ax(m_j), 'position');
                            p_col = get(h_col(m_j), 'position');
                            p_col([2, 4]) = p_ax([2, 4]);
                            set(h_col(m_j), 'position', p_col);
                            set(h_ax(m_j), 'position', p_ax);
                        end
                    end
                end
                
                % match color axes on cc&rho, dphi&dF1;
                cc_rho_ind = cellfun(@(ms) find(strcmp(ms, measures_h),1), {'cc', 'rho'}, 'un', 0);
                cc_rho_ind = [cc_rho_ind{:}];
                clims = get(h_ax(cc_rho_ind), 'clim');
                if length(cc_rho_ind) > 1
                    clims = [clims{:}];
                end                    
                if ~isempty(cc_rho_ind) && strcmp(dataToImage, 'data')
                    clims = [-1, 1]*max(abs(clims));
                    set(h_ax(cc_rho_ind), 'clim', clims);                                                                                                                
                end
            
                dphi_dF1_ind = cellfun(@(ms) find(strcmp(ms, measures_h),1), {'dphi', 'dF1'}, 'un', 0);
                dphi_dF1_ind = [dphi_dF1_ind{:}];
                clims = get(h_ax(dphi_dF1_ind), 'clim');
                if length(dphi_dF1_ind) > 1
                    clims = [clims{:}];
                end
                if ~isempty(dphi_dF1_ind) && strcmp(dataToImage, 'data')
                    nullVal = iff(strcmp(dataToImage, 'data'), 90, 0);                                        
                    clims = [-1, 1]*max(abs(clims-nullVal))+nullVal;
                    set(h_ax(dphi_dF1_ind), 'clim', clims);
                end                
                
                for bp_i = 1:length(plotsOnBottom)
                    switch plotsOnBottom{bp_i}
                        case 'NPairs',
                            Cvals = log10(numPairs{p_j,1,1});
                            cornerVals = numPairs{p_j,1,1};
                            bp_titleStr = 'log(NPairs)';
                        case 'pvals', 
                            Cvals = -log10(pvals_plot);
                            cornerVals = pvals_plot;
                            bp_titleStr = 'p-value';
                        case 'CC',
                            Cvals = log10(numPairs{p_j,1,1});
                            cornerVals = numPairs{p_j,1,1};
                            bp_titleStr = 'log(NPairs)';
                    end
                            
%                     h_wcc_N = subplot(M, N, nM_h+1,n_ac,subplotInds(nM_h+1,n_ac,nSubInds{:}));
%                     figure(1055);  h_wcc_bp(bp_i) = axes; 
                    h_wcc_bp(bp_i) = subplot(M,4, nBottomPlots_idx{bp_i}); 
                    imagesc(xs_cent, ys_cent, Cvals');
                    title(bp_titleStr); xlabel(X.labl); ylabel(Y.labl);                     
                    set(h_wcc_bp(bp_i), 'active', 'outer');
%                     putCornerNumbers(h_wcc_bp(bp_i), xs_cent, ys_cent, cornerVals, Cvals);
                    axis xy;
                    colorbar;
                        
                end
                if ~skipBigHeading
                    titleStr = [titleCase(gratingType_s) ' Gratings ' pairTypeLabels{p_j}];
%                     suptitle_2(titleStr);                            
                end
                
            end
            
            
                            
    %             text(X.vals(1),  Y.vals(1),   num2str(Wcc_NumPairs(1,1)),   'vert', 'top', 'hor', 'left', 'color', 'w');
    %             text(X.vals(1), Y.vals(end),  num2str(Wcc_NumPairs(end,1)),  'vert', 'bot', 'hor', 'left', 'color', 'k');
    %             text(X.vals(end), Y.vals(1),  num2str(Wcc_NumPairs(1,end)),  'vert', 'top', 'hor', 'right', 'color', 'k' );
    %             text(X.vals(end), Y.vals(end),num2str(Wcc_NumPairs(end,end)),  'vert', 'bot', 'hor', 'right', 'color', 'w');


%             end
            
            
        elseif doErrorBarPlots   % single 1D plot
            
%             if showPvalPlot
%                 n_ac = 2;
%                 nSubInds = {[3,1], [3,1]};
%                 pSubInds = {[3,2], [3,2]};            
%             else
%                 n_ac = 8;
%                 nSubInds = {[3,3], [3,6]};
%             end            
            line2 = '';
            doErrBar= true;
            nAx = iff(isempty(line2), 1, 2);
            
            plot2ColsIfOneLoc = false;
            
            if plot2ColsIfOneLoc && (nLoc_h == 1) && ~odd(nM_h)
%                 M = nM_h/2+1;
                M = nM_h/2;
                N = nM_h/2;
            else
%                 M = nM_h+1;
                M = nM_h;
                N = nLoc_h;                
            end
%             M = 3;
            nPairs_idx = subplotInds(M,5, [M,1], [M,3]);                
            
            cols = 'brgk';
            x_1 = xs_cent(1); x_end = xs_cent(end); 
            if ibetween(1-x_end, 0, .15), x_end = 1; end;
            if ibetween(x_1, 0, .15), x_1 = 0; end;
            
            xlims = [x_1, x_end];
            if doErrBar
                d = diff(xlims)/35;
                xlims = xlims + [-d, d];
            end
            
            h_ax = zeros(nP_h, nLoc_h, nM_h, nAx);
            h_ls = zeros(nP_h, nLoc_h, nM_h, nAx, ny);
            for p_j = 1:nP_h;                
                plot_j = 1;
                pair_str = iff(nP_h>1,  '', [pairTypes1{p_j} ' : ']);                                
                
                figure(200+p_j); clf; set(gcf, 'Name', gratingType_s);
                
                for loc_j = 1:nLoc_h
                    loc_str = iff(nLoc_h>1, '', [locations1{loc_j} ' : ']);
                    for m_j = 1:nM_h                                                                        
                        ms_str = iff(nM_h>1, '', measures_h{m_j});
                                                
                        h(plot_j) = subplot(M,N,plot_j);                        
%                         h(plot_j) = subplot(M,N,1:2);                        
                                                
                        if doErrBar
                            plot1Func = @(x, y) errorbar(x, y, data_std{p_j, loc_j, m_j});
                        else
                            plot1Func = @plot;
                        end
                                                
                        if ~isempty(line2)
                            [h_ax(p_j, loc_j, m_j,:), h_y1, h_y2] = plotyy(xs_cent, data{p_j, loc_j, m_j}, xs_cent, numPairs{p_j, loc_j, m_j}, plot1Func, @semilogy);
                        else
                            dx = diff(xlims)/150;
                            xs_plot = repmat(xs_cent, [ny,1])';
                            xs_plot = bsxfun(@plus, xs_plot, linspace(-1, 1, ny)*dx);
                            h_y1 = plot1Func(xs_plot, data{p_j, loc_j, m_j});                                
                            h_ax = gca;
                        end
                        
                        
%                         for y_i = 2:ny
%                             set(h_ax, 'nextPlot', 'add');
%                             h_y1(y_i) = plot1Func(xs_cent+dx, data{p_j, loc_j, m_j});
%                             set(h_y1(y_i), 'color', cols(y_i))
%                         end
                            
%                         [h_ax, h_y1, h_y2] = plotyy(xs_cent, data{p_j, loc_j, m_j}, xs, -log10(pvals{p_j, loc_j, m_j}), eb, @plot);

                        % x limits
                        set(h_ax, 'xlim', xlims);

                        % y1 limits/ticks                        
                        ydata = data{p_j, loc_j, m_j};
                        if doErrBar
                            ydata = [ydata + data_std{p_j, loc_j, m_j}; ydata - data_std{p_j, loc_j, m_j}];
                        end
                        switch measures_h{m_j}
                            case {'dot', 'cc', 'rho'}, 
                                nullMean = 0;
                                ylims = [roundToNearest(min([nullMean; ydata(:)]), .025, 'down'), ...
                                         roundToNearest(max( ydata(:)    ), .025, 'up') ];                                
                                yticks_lims = [roundToNearest(ylims(1), .1, 'down'), ...
                                               roundToNearest(ylims(2), .1, 'up')];
                                yticks = yticks_lims(1) : .1 : yticks_lims(end);
                                fmt = '%.1f';
                            case {'dphi', 'dF1'}, 
                                nullMean = 90;
                                ylims = [roundToNearest(min(ydata(:)), 2.5,       'down'), ...
                                         roundToNearest(max([nullMean; ydata(:)]), 2.5, 'up')];                                            
                                yticks_lims = [roundToNearest(ylims(1), 10, 'down'), ...
                                               roundToNearest(ylims(2), 10, 'up')];
                                d_tck = iff(diff(yticks_lims) <80, 10, 20);
                                yticks = yticks_lims(1) : d_tck : yticks_lims(end);
                                fmt = '%d';
                        end                        
%                         set(h_ax(1), 'ylim', ylims, 'ytick', yticks, 'ytickLabel', num2str(yticks', fmt) );
                        set(h_ax(1), 'ylim', ylims );
                        drawHorizontalLine(nullMean, 'linestyle', ':');
                        set(h_y1, 'marker', 'o', 'markersize', 3);
                        if addAxesLabelsOnEachPlot
                            set(get(h_ax(1),'Ylabel'),'String',measureLabels_h{m_j}) ;
                        end
                        
                        % y2 limits/ticks                        
                        if ~isempty(line2)
                            yticks2 = [floor(log10(min(numPairs{p_j, loc_j, m_j})) ), ceil( log10(max(numPairs{p_j, loc_j, m_j})) )];
                            set(h_ax(2), 'ylim', 10.^[yticks2], 'ytick', 10.^[yticks2(1):yticks2(2)]);
                            set(h_y2, 'marker', '.');
                            set(get(h_ax(2),'Ylabel'),'String','N');
                        end
                        set(h(plot_j), 'active', 'outer');
                                                
%                         title([pair_str ' ' loc_str, ms_str ]);

%                         if m_j > 2
                            xlabel(X.labl);
%                         end                        
                        
                        
                        pairsOfCols = nchoosek(1:ny,2);                        
                        nPrs = size(pairsOfCols,1);
                        sig_Ps = ones(nx, nPrs);
                        ms_vals = cellfun(@(v) v{m_j}, setVals(p_j,:,:), 'un', 0);
                        ms_vals = reshape(ms_vals, [nx, nPrs]);
                        ms_idx = cellfun(@(v) v{m_j}, setIdx(p_j,:,:), 'un', 0);
                        ms_idx = reshape(ms_idx, [nx, nPrs]);
                        for xi = 1:nx
                            for pr_i = 1:nPrs
                                idx = pairsOfCols(pr_i,:) ;
                                [v1, v2] = deal(ms_vals{xi,idx});
                                if length(v1) > 1 && length(v2) > 1
                                    sig_Ps(xi,pr_i) = ranksum(v1, v2);
                                    [tmp, sig_Ps(xi,pr_i)] = ttest2(v1, v2);                                
                                end
                                    
                            end
                        end               
                        if showSigTextOnAxes
                            sig_Ps = min(sig_Ps, [], 2);
                            if any(sig_Ps(:) < .05)
                                sig_pos = find(sig_Ps(:) < .05);
                                for jj = 1:length(sig_pos)
                                    L = num2str(-log10(sig_Ps(sig_pos(jj))), '%.1f');
    %                                 {'\fontsize{15}*', '\fontsize{9}' L
    %                                 text(xs_cent(sig_pos(jj)), ylims(1)+diff(ylims)*.9, *, 'fontsize', 20, 'hor', 'cent')

                                    text(xs_cent(sig_pos(jj)), ylims(1)+diff(ylims)*.99, L, 'fontsize', 6, 'hor', 'cent', 'vert', 'top')
                                    
                                end                            
                            else
    %                             annotation('textbox', [.03, .08, .1, .1], 'String', 'X');
                                text(xlims(1) + diff(xlims)*.03 , ylims(1)+diff(ylims)*.08, 'X', 'fontsize', 10, 'hor', 'cent');
                            end
                        end
                    
                        
                        plot_j = plot_j+1;
                    end
                end
                
                cc_rho_ind = cellfun(@(ms) find(strcmp(ms, measures_h),1), {'cc', 'rho'}, 'un', 0);
                cc_rho_ind = [cc_rho_ind{:}];
                if length(cc_rho_ind) > 1
                    matchAxes('Y', h(cc_rho_ind));
                end
                dphi_dF1_ind = cellfun(@(ms) find(strcmp(ms, measures_h),1), {'dphi', 'dF1'}, 'un', 0);
                dphi_dF1_ind = [dphi_dF1_ind{:}];
                if length(dphi_dF1_ind) > 1
                    matchAxes('Y', h(dphi_dF1_ind));
                end                
                
%                 strcmp(
                
                figure(1056); clf;
%                 h_N = subplot(M,5, nPairs_idx);
                h_N = subplot(1,5, [1:3]);
                semilogy(xs_cent, numPairs{p_j,1,1}, '.-');
                allNs = numPairs{p_j,1,1}(:);
                ylims = [max(min(allNs),1), max(allNs)];
                yticks = 10.^[floor(log10(ylims(1))): ceil(log10(ylims(2)))];
                set(gca, 'ytick', yticks)
                ylim([yticks(1), yticks(end)])
                
                ylabel(['Number of Pairs']);
                xlabel(X.labl); 
                set(h_N, 'xlim', xlims);
                set(h_N, 'active', 'outer');    
%                 legend( legendarray( [Y.labl ' = '], Y.vals), 'location', 'best' )
                p = get(gca, 'position');
                if ~iscell(Y.vals) 
                    legend( legendarray( [Y.labl ' = '], Y.vals), 'location', 'bestoutside' )
                else % special case
                    leglabs = arrayfun(@(pre, post) sprintf('%.2f < %s < %.2f', pre, Y.labl, post), Y.vals{1}, Y.vals{2}, 'un', 0);
                    legend(leglabs, 'location', 'bestoutside' )                    
                end
                set(gca, 'position', p)
%             if showPvalPlot
%                 h_wcc_P = subplot(3,n_ac,subplotInds(nM_h+1,n_ac,pSubInds{:})); imagesc(xs_cent, ys, -log10(pvals_plot)');
%                 title([p_s '-log(p)']); xlabel(X.labl); ylabel(Y.labl); colorbar;
%                 putCornerNumbers(h_wcc_P, X, Y, pvals_plot, -log10(pvals_plot));
%                 set(h_wcc_P, 'active', 'outer');                
%             end            
                titleStr = [titleCase(gratingType_s) ' Gratings ' pairTypeLabels{p_j}];                
%                 suptitle_2(titleStr);   

            end
         
            
        end
        
        if saveDataToFile
            paramSearchFilename = [CatV1Path gratingType_s 'ParameterSearch1_' cmpType_s '_' pairTypes_str '.mat'];
            save(paramSearchFilename, 'X', 'Y', 'Z', 'measures_h', 'data', 'data_std', 'numPairs', 'pvals', 'pvalLabels_h');
        end
        
        keyboard;
    end    
    
    
    if makeBarGraphOfResults 
        if length(pairTypes) < 2
            error('Need to have both Wcc and Bcc active');
        end
%         rescaleDphiDistribution = false;
        rescaleBccDistribution = false;
            
        % calculate:
        [toDisplay,nResamples, colorCoding] = deal(toDisplayCur, 1, '[none]');         
        basicInputs = {toDisplay, measures{1}, measures{2}, pairTypesActive, nResamples, locationsActive, measuresActive, colorCoding};        
        filterInitVals = cellfun(@(C) C{3}, filterArgs, 'un', 0);        
        filterInputs_all = cell2struct(filterInitVals, filterArgNames, 2);
        
        % initialize: set values        
%         filterInputs = enableFilter(filterInputs0_nofilter, 'min_rep_cc_p');
%         filterInputs.minRepPval_cc_p = 2;

        % after initial compilation  set X & Y filters on.                        
                
        filterInputs_allFrac = enableFilter(filterInputs_all, 'minFracR_lo');        
        filterInputs_allFrac.minFracR_lo = .65; %iff(gratingType == 1, .65, .65);
        
        filterInputs_CC = enableFilter(filterInputs_allFrac, 'maxF1oDC_cmp_hi');
        filterInputs_CC.maxF1oDC_cmp_hi = .7;
        
        filterInputs_SS = enableFilter(filterInputs_allFrac, 'minF1oDC_cmp_lo');        
        filterInputs_SS.minF1oDC_cmp_lo = .75;        
        filterInputs_SS = enableFilter(filterInputs_SS, 'maxF1oDC_cmp_hi');        
        filterInputs_SS.maxF1oDC_cmp_hi = iff(gratingType == 1, 1.7, 1.5);

        
        
        filterInputs_all_C = struct2cell(filterInputs_all);
%         filterInputs_allFrac_C = struct2cell(filterInputs_allFrac);
        filterInputs_CC_C = struct2cell(filterInputs_CC);
        filterInputs_SS_C = struct2cell(filterInputs_SS);
        
        
        [distMean_all, distStd_all, distStderr_all, distN_all, distP_all, binN_all, tmpBinDN_all, sVals_all, sIdx0_all] = updateFilters(basicInputs{:}, filterInputs_all_C{:});
%         [distMean_frac, distStd_frac, distStderr_frac, distN_frac, distP_frac, binN_frac, tmpBinDN_frac, sVals_frac, sIdx0_frac] = updateFilters(basicInputs{:}, filterInputs_allFrac_C{:});
        [distMean_CC, distStd_CC, distStderr_CC, distN_CC, distP_CC, binN_CC, tmpBinDN_CC, sVals_CC, sIdx_CC] = updateFilters(basicInputs{:}, filterInputs_CC_C{:});
        [distMean_SS, distStd_SS, distStderr_SS, distN_SS, distP_SS, binN_SS, tmpBinDN_SS, sVals_SS, sIdx_SS] = updateFilters(basicInputs{:}, filterInputs_SS_C{:});

        x_bin = 1:4;
                
        
        figure(105); clf; set(gcf, 'Name', gratingType_s);

        for m_j = 1:nM_h
            h(m_j) = subplot(4,1,m_j);
            
            vals = [sVals_all(2, 1, m_j), sVals_all(1, 1, m_j), sVals_CC(1, 1, m_j), sVals_SS(1, 1, m_j)];
            Y = [distMean_all(2, 1, m_j), distMean_all(1, 1, m_j), distMean_CC(1, 1, m_j), distMean_SS(1, 1, m_j)];
            E = [distStderr_all(2, 1, m_j), distStderr_all(1, 1, m_j), distStderr_CC(1, 1, m_j), distStderr_SS(1, 1, m_j)];
            Pv = [distP_all(2, 1, m_j), distP_all(1, 1, m_j), distP_CC(1, 1, m_j), distP_SS(1, 1, m_j)];
            
            col = [.8 .8 .8];
            baseValue = mMullMeans(m_j);            
            hbar = bar(x_bin, Y, 'facecolor', col, 'baseValue', baseValue); hold on;
            baseline_handle = get(hbar,'BaseLine');
            set(baseline_handle,'LineStyle',':')
            set(gca, 'xtick', [])
            
            idx_toTest = 2:x_bin(end);
            valsCtrl = vals{1};
            pvals = zeros(1,length(x_bin));            
            for pv_i = idx_toTest;
                if strcmp(measures_h{m_j}, 'dphi')
                    p_chsqr_id = find( strcmp('p_{X}', curPvalLabels{1,1,m_j}), 1 );            
                    pvals(pv_i) = Pv{pv_i}(p_chsqr_id);
                else
                    pvals(pv_i) = ranksum(valsCtrl, vals{pv_i});
                end
                
            end
            
%                     idx_dphi = find(find(measuresActive)==3, 1);
%                     p_chsqr_id = find( strcmp('p_{X}', curPvalLabels{1,1,idx_dphi}), 1 );            
            
            
            errorb(x_bin, Y, E, 'linewidth', 1.5, 'barwidth', 1, 'color', 'k')
            if any(strcmp(measures_h{m_j}, {'cc', 'rho'}))
                set(gca, 'ytick', [0 :.1:.3])
                if gratingType == 1
                    ylim([-.05, .38]);                                    
                else
                    ylim([-.03, .32]);
                end                     
                sgn = 1;

            else                
                if gratingType == 1
                    set(gca, 'ytick', [60, 70, 80, 90])
                    ylim([55, 95]);
                else
                    set(gca, 'ytick', [60, 70, 80, 90])
                    ylim([60, 93]);
                end
                sgn = -1;
            end
            xlim([.5, x_bin(end)+.5]);                
%             pairings = [1,2; 2,3; 1,3];
%             levels = [1, 1, 2];
%             for pr_i = 1:size(pairings,1),
%                 pval_pairs(pr_i) = ranksum(vals{pairings(pr_i,1)}, vals{pairings(pr_i,2)});
%             end
%             drawPvalLines(x_bin(pairings), Y(pairings)+sgn*E(pairings), pval_pairs, levels, sgn, false)
            
            drawPvalStars(x_bin(idx_toTest), Y(idx_toTest)+sgn*E(idx_toTest), pvals(idx_toTest), sgn, true)
            
            if (m_j == nM_h)
                xtickLabels = {'Ctrl', 'All', 'C/C', 'S/S'};
%                 xtickLabels = arrayfun(@(i) ['(' char(i+64) ')'], x_bin, 'un', 0);
                set(gca, 'xtick', x_bin, 'xtickLabel', xtickLabels)
            end
    
            %                         title([p_s ' ' locations1{loc_j} ' : ' measures_h{m_j} ]);
            
%             title(measureLabels_h{m_j}) ;
%             plot_j = plot_j+1;
        end
        
        3;
    end
    
end


function drawPvalLines(xval_pairs, yval_pairs, pval_pairs, levels, sgn, justStars)
    assert(any(sgn == [-1, 1]));

    nPairs = size(xval_pairs,1);
    ylims = get(gca, 'ylim');
    gap_frac = .1;
    if sgn == 1
        y_start = max(yval_pairs(:))*(1 + gap_frac);
        y_end   = ylims(2) - diff(ylims)*.05;
    else
        y_start = min(yval_pairs(:));
        y_start = y_start - (90-y_start)*gap_frac;
        y_end   = ylims(1) + diff(ylims)*.05;
    end    
    
    hTot = abs(y_end - y_start); %hTot = diff(ylims)*(.1);
    nLevels = length(unique(levels));
    hLev = hTot/nLevels;
    txtRatio = .7;
    hBot = y_start + sgn*(0:nLevels-1)*hLev;
    
    dx = min(diff(unique(xval_pairs(:))))/20;
    
    for pi = 1:nPairs
        x1 = xval_pairs(pi,1);
        x2 = xval_pairs(pi,2);
        
        x1 = x1+dx;
        x2 = x2-dx;
        
        xT = (x1+x2)/2;
        y1 = hBot(levels(pi));
        y2 = y1+sgn*(1-txtRatio)*hLev;
%         yT = y1+(txtRatio)*hLev
%         y1 = yval_pairs(pi,1);
%         y2 = yval_pairs(pi,2);
        p  = pval_pairs(pi);
        
        X = [x1, x1, x2;
             x1, x2, x2];
        Y = [y1, y2, y2;
             y2, y2, y1];
        pval = pval_pairs(pi);
        if justStars
            str = pval2stars(pval);
        else
            fmt = iff(pval < .01, '%.1e', '%1.2g');
            str = sprintf(['p = ' fmt], pval);
            str = strrep(str, 'e-0', 'e-');
            str = strrep(str, 'e-0', 'e-');
        end
        
        hL(pi,:) = line(X, Y, 'color', 'k');     
        vPos = iff(sgn == 1, 'bottom', 'top');
%         vPos = iff(sgn == 1, 'baseline', 'cap');
        
        hText(pi) = text(xT, y2, str, 'fontsize', 10, 'horizontal', 'center', 'vertical', vPos);        
    end


end


function drawPvalStars(xvals, yvals, pvals, sgn, justStars)
    assert(any(sgn == [-1, 1]));

    nPairs = length(xvals);
    ylims = get(gca, 'ylim');
    gap_frac = .1;
    if sgn == 1
        y_start = max(yvals(:))*(1 + gap_frac);
        y_end   = ylims(2) - diff(ylims)*.05;
    else
        y_start = min(yvals(:));
        y_start = y_start - (90-y_start)*gap_frac;
        y_end   = ylims(1) + diff(ylims)*.05;
    end    
    
    hTot = abs(y_end - y_start); %hTot = diff(ylims)*(.1);        
    
    for pi = 1:nPairs
        x = xvals(pi);
%         if (sgn == -1)
            dy = -diff(ylims)*.02;
%         elseif (sgn == 1)
%             dy = -diff(ylims)*.02;
%         end
        y = yvals(pi) + dy;
        pval = pvals(pi);
        
        if justStars
            str = pval2stars(pval);
        else
            fmt = iff(pval < .01, '%.1e', '%1.2g');
            str = sprintf(['p = ' fmt], pval);
            str = strrep(str, 'e-0', 'e-');
            str = strrep(str, 'e-0', 'e-');
        end        
        vPos = iff(sgn == 1, 'bottom', 'top');
%         vPos = iff(sgn == 1, 'baseline', 'cap');
        
        hText(pi) = text(x, y, str, 'fontsize', 10, 'horizontal', 'center', 'vertical', vPos);
    end

end


function s = pval2stars(pval)
    pvalInts  = [ .001,  .01,  .05,    1 ];
    pvalStrs  = {'***',  '**', '*' , 'ns'};

    id = find(pval <= pvalInts, 1, 'first');
    s = pvalStrs{id};
    
%     if (pval > pvalInts(1))
%         s = 'ns';
%     elseif (pval <= pvalInts(1)) && (pval > pvalInts(2))
%         s = '*';
%     elseif (pval <= pvalInts() && (pval > .005)
%         s = '**';
%     elseif (pval <= .001) %&& (pval > .001)
%         s = '***';
%     elseif (pval <= .001)
%         s = '****';
%     end
       
end


function  putCornerNumbers(ax, xs, ys, n_txt, n_col)
    axes(ax)
%     decfmt = '%.1e';    
    mn = min(n_col(:)); mx = max(n_col(:));
    scl = (n_col-mn)/(mx-mn);
    bw = @(ni) iff(ibetween(ni, .35, .68), 'k', 'w');
    
    text(xs(1),   ys(1),   n2str(n_txt(1,1)),   'vert', 'top', 'hor', 'left', 'color', bw(scl(1,1)));
    text(xs(1),   ys(end), n2str(n_txt(end,1)),  'vert', 'bot', 'hor', 'left', 'color', bw(scl(end,1)));
    text(xs(end), ys(1),   n2str(n_txt(1,end)),  'vert', 'top', 'hor', 'right', 'color', bw(scl(1,end)) );
    text(xs(end), ys(end), n2str(n_txt(end,end)),  'vert', 'bot', 'hor', 'right', 'color', bw(scl(end,end)));
end

function X = adjustImageNans(X)
    f = 0.03;
    mn = min(X(:));
    mx = max(X(:));
    nan_val = mn-f*(mx-mn);
    X(isnan(X)) = nan_val;    
end

% function getSamplePvalueFromData(x)
%     n = length(x);    
%     x_pn = [x(:); -x(:)];
% 
%     
%     
%     bootstrp(1000, @(x) 
% 
% 
% 
% 
% end

function s = n2str(n)    
    if n == round(n)
        s = num2str(n, '%d');
    else
        if n >= .1
            s = num2str(n, '%.1f');
        else
            s = num2str(n, '%.0e');
            
            [s_i] = regexp(s, 'e-0*', 'once');
            if ~isempty(s_i)
                s = regexprep(s, 'e-0*', 'e-');
            end
            
        end
    end
%     fmt = iff(any( n(:) ~= round(n(:))), decfmt, '%d');
end

function filterInputs = enableFilter(filterInputs, thName, isOn)
    if nargin < 3
        isOn = true;
    end
    flds = fieldnames(filterInputs);    
    filterHeadingIdxs = find( cellfun(@(s) ~isempty(strfind(s, 'Filter')), flds) );
    
    if ischar(thName)
        thName = {thName};
    end
    for th_i = 1:length(thName);
        th_idx = find(strcmp(thName{th_i}, fieldnames(filterInputs) ));
        if isempty(th_idx)
            error('Invalid fieldname');
        end
        curHeaderIdx = filterHeadingIdxs( find(filterHeadingIdxs < th_idx, 1, 'last') );
        offset = th_idx - curHeaderIdx;
        filterInputs.(flds{curHeaderIdx})(offset) = logical(isOn);    
    end        
end
                                

function filterInputs = setFilterValue(filterInputs, varName, varVals, varValIdx)
    if iscell(varName)
        for vi = 1:length(varName)
            filterInputs.(varName{vi}) = element(varVals{vi}, varValIdx);
        end
    else
        filterInputs.(varName) = element(varVals, varValIdx);
    end
end


%     showPairOptions = repmat({tf}, 1, length(pairTypesAvailable)); pairTF(~pairTypesAvailable) = {false}; % make these automatically not-shown    
%     locTF   = repmat({tf}, 1, length(locationsAvailable)); locTF(~locationsAvailable)  = {false}; % make these automatically not-shown        
%     msTF   = repmat({tf}, 1, length(measuresAvailable));  msTF(~measuresAvailable)    = {false}; % make these automatically not-shown    
%     pairInit = true(1, length(pairTypesAvailable)); pairInit(~pairTypesAvailable) = false; % make these automatically not-shown    
%     locInit  = true(1, length(locationsAvailable)); locInit(~locationsAvailable)  = false; % make these automatically not-shown    
%     msInit   = true(1, length(measuresAvailable));  msInit(~measuresAvailable)    = false; % make these automatically not-shown    


%{
%                 p_bcc = find(strcmp(pairTypes, 'Bcc'));
%                 p_wcc = find(strcmp(pairTypes, 'Wcc'));
% 
%                 wcc_id = single(color(pairIdxs{p_wcc}));
%                 if any(strcmp(measures{ms}, nPh_sensitiveMeasures))
%                     wcc_id = [wcc_id, [pairData(pairIdxs{p_wcc}).n_phases]']; %#ok<AGROW>
%                 end
% %                 wcc_nph = [pairData(setIdx{p_wcc,loc,ms}).n_phases];
%                 [wcc_uph, wcc_ph_count] = uniqueCount(wcc_id, 'rows');
% 
%                 bcc_id = single(catIds(pairIdxs{p_bcc}));
% %                 bcc_id = single(colorMask(pairIdxs{p_bcc}));
% %                 bcc_id = [pairData(pairIdxs{p_bcc}).SCtype]';
%                 if any(strcmp(measures{ms}, nPh_sensitiveMeasures))
%                     bcc_id = [bcc_id, [pairData(pairIdxs{p_bcc}).n_phases]']; %#ok<AGROW>
%                 end
%                 [bcc_uph, bcc_ph_lists] = uniqueList(bcc_id, 'rows');
%                 bcc_ph_count = cellfun(@length, bcc_ph_lists);
% %                 rand('state', 0);
%                 if isempty(setxor(wcc_uph, bcc_uph, 'rows'))
%                     [bcc_matched_count, eff] = matchNumbersToRatios(bcc_ph_count, wcc_ph_count/sum(wcc_ph_count));
%                     bcc_nToRemove = bcc_ph_count - bcc_matched_count;
%                     for i = 1:length(bcc_ph_lists)
%                         idx = randperm( bcc_ph_count(i)); idx = idx(1:bcc_nToRemove(i));
%                         bcc_ph_lists{i} = bcc_ph_lists{i}(idx(:));
%                     end
%                     remove_idx = cat(1, bcc_ph_lists{:});
%                     mask(pairIdxs{p_bcc}(remove_idx)) = false;
%                 end
%             end
%}

%{
% Steps to add a new limit (before automated!!)
% 1) Beginning of file: Add to pairFilterVars/locFilterVars/valFilterVars
% 
% 2) In function updateFilters:
% 	- add to varargin:
% 	- add to filterActive (in SAME ORDER as was in beginning of file)
% 	- add to th
% 
% 3) At end, in control panel vars:
% 	Groups: add to limitByVars
% 	Add to args list
%   
% 4) for heat maps: add to list of inputs for updateFilters call
%
% 5) make sure to add the necessary fields in generateGratingComparisonDatafiles.m !



old (un-automated) generation of filter args for "manipulate"




    pairFilters1_labels  = {'rep_p_av', 'n_phases'};    nCf = length(pairFilters1_labels);
    pairFilters1_depVars = { {{}, {'min_rep_cc_p'}}, {{}, {'n_phases'}} };
    
    pairFilters2_labels = {'S/C', 'F1/DC_pref_lo', 'F1/DC_pref_hi'};   nPf1 = length(pairFilters2_labels);
    pairFilters2_depVars ={ {{}, {'SCtype'}}, {{}, {'minF1oDC_pref_lo'}}, {{}, {'minF1oDC_pref_hi'}} };
    
    pairFilters3_labels = {'Animal', 'Penetration', 'Location'};   nPf2 = length(pairFilters3_labels);        
    pairFilters3_depVars = { {{}, {'animalCmp'}} , {{}, {'penetrCmp'}}, {{}, {'locCmp'} } };
            
    locFilters1_labels = {'F1/DC_cmp_lo', 'F1/DC_cmp_hi', 'sumPhs'};  nLf1 = length(locFilters1_labels);
    locFilters1_depVars ={ {{}, {'minF1oDC_cmp_lo'}}, {{}, {'minF1oDC_cmp_hi'}}, {{}, {'minSumPhs'}} };

    locFilters2_labels = {'minFracR_lo', 'minFracR_hi'};  nLf2 = length(locFilters2_labels);
    locFilters2_depVars ={ {{}, {'minFracR_lo'}}, {{}, {'minFracR_hi'}} };
    
    setFilters_labels  = { {'pval'} };   nSf = length(setFilters_labels);
    setFilters_depVars = { {}, {'pval'} } ;
    

    filterArgs = {...
             {'pairFilter1', repmat(tf, nCf,1), true(nCf,1), pairFilters1_depVars, pairFilters1_labels},  ...
                  {'min_rep_cc_p', logPvalRng, th.min_rep_cc_p}, ...                  
                  {'n_phases', allNPhases_S, allNPhases_S{1}}, ...
             {'pairFilter2', repmat(tf, nPf1,1), true(nPf1,1), pairFilters2_depVars, pairFilters2_labels}, ...
                 {'SCtype', SCtypesStr, SCtypesStr{1}}, ...
                 {'minF1oDC_pref_lo', rng01, th.minF1oDC_pref_lo}, ...
                 {'minF1oDC_pref_hi', rng01, th.minF1oDC_pref_hi}, ...
             {'pairFilter3', repmat(tf, nPf2,1), true(nPf2,1), pairFilters3_depVars, pairFilters3_labels}, ...
                 {'animalCmp', locCmpStrs, locCmpStrs{1}}, ...
                 {'penetrCmp', locCmpStrs, locCmpStrs{1}}, ...
                 {'locCmp', locCmpStrs, locCmpStrs{1}}, ...
             {'locFilter1', repmat(tf, nLf1,1), true(nLf1,1), locFilters1_depVars, locFilters1_labels}, ...
                 {'minF1oDC_cmp_lo', rng01, th.minF1oDC_cmp_lo}, ...
                 {'minF1oDC_cmp_hi', rng01, th.minF1oDC_cmp_hi}, ...
                 {'minSumPhs', [0:.25:20], th.minSumPhs}, ...
             {'locFilter2', repmat(tf, nLf2,1), true(nLf2,1), locFilters2_depVars, locFilters2_labels}, ...
                 {'minFracR_lo', rng01, th.minFracR_lo}, ...
                 {'minFracR_hi', rng01, th.minFracR_hi}, ...
             {'valFilter1', repmat(tf, nSf,1), true(nSf,1), setFilters_depVars, setFilters_labels}, ...
                {'pval',  pvalsRng, th.pval} ...
            };

%}


%{
delta phi with first/2nd /... harmonics
discrete dF1 vs delta phi.
discrete dF1+dF2 vs delta phi
see how many harmonics it takes for them to agree

-- which smoothing functions work best
	all harmonics/spline/etc.

distribution of power vs harmonics

std/mean, # harmonics, smoothed-non-smoothed
%}



%{

        if multipleXth
            filterInputs = enableFilter(filterInputs, X.thName{1});
            filterInputs = enableFilter(filterInputs, X.thName{2});
        else
            filterInputs = enableFilter(filterInputs, X.thName);
        end
        if multipleYth
            filterInputs = enableFilter(filterInputs, Y.thName{1});
            filterInputs = enableFilter(filterInputs, Y.thName{2});
        else
            filterInputs = enableFilter(filterInputs, Y.thName);
        end

        if multipleXth
            filterInputs.(X.thName{1}) = element(xs{1}, x_i);
            filterInputs.(X.thName{2}) = element(xs{2}, x_i);
        else
            filterInputs.(X.thName) = element(xs, x_i);
        end

        if multipleYvals,
            if multipleYth                        
                filterInputs.(Y.thName{1}) = element(ys{1}, y_i);
                filterInputs.(Y.thName{2}) = element(ys{2}, y_i);
            else
                filterInputs.(Y.thName) = element(ys, y_i);
            end
        end
        if ~isempty(Z)                        
            if multipleZth
                filterInputs.(Z.thName{1}) = element(zs{1}, z_i);
                filterInputs.(Z.thName{2}) = element(zs{2}, z_i);
            else
                filterInputs.(Z.thName) = element(zs, z_i);
            end
        end

%}

%{
%%
global  minFrac_Bcc minFrac_Wcc
minFrac_Bcc = pairData.loc_minFracOfMaxes(idxMtx(Bcc_pairIdxs));
minFrac_Wcc = pairData.loc_minFracOfMaxes(idxMtx(Wcc_pairIdxs));
h = hist2({minFrac_Wcc, minFrac_Bcc}, [0:.0333:1], 'norm', 'stairs')
xlim([0 1])
xlabel('maxMinFracR'); ylabel('P(maxMinFracR)')
legend({'Within-site', 'Between-site'}, 'location', 'NW')
%}

%{
% find a 3 cells: 2 from the same group that have high minFracR, and of the sam
%}
