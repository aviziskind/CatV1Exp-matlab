
%%
redoAllCorticalStates = 1;
if redoAllCorticalStates
    %%
    [allGids,allCellIds] = getAllGids;
%     [allGids,allCellIds] = getAllGids('f');
    nCells = length(allGids);
    %%
    S_or = load('cellsGroups_GLFcuw8_grating_dOr.mat');
    oG = [S_or.gratingGroups_dOr.Gid];
    oOr = [S_or.gratingGroups_dOr.stimOrdered];

    S_sf = load('cellsGroups_GLFcuw8_grating_dSf.mat');
    sG = [S_sf.gratingGroups_dSf.Gid];
    sOr = [S_sf.gratingGroups_dSf.stimOrdered];

    Gid_skip = [oG(oOr == 1), sG(sOr == 1)];
    
      %%  
    i_start = 1;
    fprintf('STARTING...\n');
    for i = i_start:nCells
%         if any(allGids(i) == Gid_skip)
%             continue;
%         end
       fprintf('%d / %d. Gid = %d. cellId = %d...', i, nCells, allGids(i), allCellIds(i));
       fprintf('sin... '); [fit_func, fitS_C{i}] = getCorticalState(allGids(i), allCellIds(i), 'sin', 8);
       fitS_C{i}.Gid = allGids(i);
       fitS_C{i}.cellId = allCellIds(i);
    %    fprintf('gauss... '); getCorticalState(allGids(i), allCellIds(i), 'gauss', 8);
       fprintf('\n');
    end
end

%%
S_fgCells = load( getFileName('indiv', 'movie_fg') ); 
% S_ori_cells = load( getFileName('indiv', 'grating_dOr') );
% S_spf_cells = load( getFileName('indiv', 'grating_dSf') );
% S_indivCells = mergeStructs(S_fgCells, S_ori_cells, S_spf_cells);
S_indivCells = S_fgCells;
% [allGids, allCellIds] = getAllGids('f');

%%
cells_filename = getFileName('osps');
S_cells = load(cells_filename);
allGids = [S_cells.allCells.Gid];
allCellIds = [S_cells.allCells.cellId];

%%
R_corr_all = {S_cells.allCells.R_corr};

allClasses = cellfun(@class, R_corr_all, 'un', 0);

idx_bad = find(strcmp(allClasses, 'MException'));
idx_ok = find(~strcmp(allClasses, 'MException'));

badGids = [S_cells.allCells(idx_bad).Gid];
badCellIds = [S_cells.allCells(idx_bad).cellId];
SGids = [S_cells.allCells.Gid];
ScellIds = [S_cells.allCells.cellId];

%%
R_all = [R_corr_all{idx_ok}];
flds = fieldnames(R_all(1));
for i = 1:length(flds)
    grpS.(flds{i}) = [R_all.(flds{i})];
end
%%
Gids_ok = allGids(idx_ok);
cellIds_ok = allCellIds(idx_ok);

%%
for i = 1:length(R_corr_all)
    
    
end

%%
% [allGids_o, allCellIds_o] = getAllGids('o');
[allGids_f, allCellIds_f] = getAllGids('f');
[allGids_s, allCellIds_s] = getAllGids('s');
addMUs = 0;
mu_cellId = 100;
onlyMUs = 1;
if addMUs
    uGids_f = unique(allGids_f);
    allGids_f = [allGids_f(:); uGids_f(:)];
    allCellIds_f = [allCellIds_f(:); ones(size(uGids_f))*mu_cellId ];

    uGids_s = unique(allGids_s);
    allGids_s = [allGids_s(:); uGids_s(:)];
    allCellIds_s = [allCellIds_s(:); ones(size(uGids_s))*mu_cellId ];
elseif onlyMUs
    uGids_f = unique(allGids_f);
    allGids_f = [uGids_f(:)];
    allCellIds_f = [ones(size(uGids_f))*mu_cellId];

    uGids_s = unique(allGids_s);
    allGids_s = [uGids_s(:)];
    allCellIds_s = [ones(size(uGids_s))*mu_cellId];
    
end
allGids    = [allGids_f(:); allGids_s(:)];
allCellIds = [allCellIds_f(:); allCellIds_s(:)];


%%
% idx_err = find(strncmp(allClass, 'ME', 2))
% idx_corr = find( [grpS.p_Rcorr_av_gain] < .05  & ScellIds > 0 );
% idx_do = idx_corr;

% Gids_do = SGids(idx_do);
% cellIds_do = ScellIds(idx_do);
%%
%%
curResponseType('gainCorrected');

%%
allR_corr_S= {};
allCCs = []; allCCs_f = [];
% for i = 1:length(idx_do)
progressBar('init-', length(allCellIds), 30)
%%
idx_f = flashedOrDrifting( allGids ) == 1 & allCellIds > 0;
idx_d = flashedOrDrifting( allGids ) == 2 & allCellIds > 0;
% allGids_f = allGids(idx_f);
% allGids_f = allGids(idx_f);
%%
% idx = find(p_Rcorr_av_gain < .0001);
% allR;
[allGids, allCellIds] = deal(allGids_s, allCellIds_s);
%%
uGids = unique(allGids);
for i = 1:length(uGids)
    fprintf('%d - %d\n', i, uGids(i))
    dbGetCellSpkStimHists(uGids(i), 0);
end
   
%%  
for i = 1:length(allCellIds)
% for i = 1:length(idx)
    %%
    
%     i = idx_err(2);  
%     Gid = allGids(i); cellId = allCellIds(i);
%     Gid = 2288; cellId = 3;
% i = i+1;
%     Gid = allGids(idx(i)); cellId = allCellIds(idx(i));
    Gid = allGids((i)); cellId = allCellIds((i));

    
    if any(Gid == Gid_skip)
        fprintf('(%d) Skipping Gid = %d\n', i, Gid)
        continue;
    end
    if any(cellId == 0)
        fprintf('(%d) Skipping multiunit\n', i);
    end
%     fprintf('\n=========================\n (%d/%d). Gid = %d, cellId = %d\n', i, length(allCellIds), Gid, cellId)
%     fitS
%
%     frmLength = getFrameLength('Gid', Gid);
%     redo = 1; %frmLength == 100;
%     Gid = Gids_ok(idx(1));  cellId = cellIds_ok(idx(1));
    
    %%
    if flashedOrDrifting(Gid) == 1
        [PSTHdata, stats] = getPSTHforCell(Gid, cellId);
        LR_bins = PSTHdata.timeWindow_bins;
        [L_bin_best, R_bin_best, windowProfile_best] = deal(LR_bins(1), LR_bins(2), PSTHdata.windowProfile);
% 
%         varname = getName('celldata', Gid, cellId);
%         celldata = S_indivCells.(varname);
%         PSTH = celldata.PSTH;
% 
%         [L_bin_best, R_bin_best, windowProfile_best] = getLRbins_windowprofile('best', PSTH);
    else
        [L_bin_best, R_bin_best, windowProfile_best] = deal(1,1,[]);
%         
    end
    %%
    fprintf('\n========================= (%d/%d), Gid = %d, cellId = %d [%d, %d]\n', i, length(allCellIds), Gid, cellId, L_bin_best, R_bin_best)

    %     Merr = [];
%     try

        getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_ph', 'osp_full'}, 1);
        getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_ph_oe', 'osp_ph_hoe'}, 1);
        3;
        %%
%         R_1 = decompress( Si.R );
%         Rf_1 = osp_full;
%         
%         [R_2, Rf_2, stats2] = getResponseCorrectedForCorticalState(Gid,cellId, [], 'all', Si.corticalGainParams);
%         assert(isequal(Si.corticalGainParams, stats2.corticalGainParams));
%         [cc_osp] = corr(R_1(:), R_2(:));
%         [cc_ospf] = corr(Rf_1(:), Rf_2(:));
%         allCCs(i) = cc_osp;
%         allCCs_f(i) = cc_ospf;
%         %%
%         figure(88); clf; subplot(2,1,1); plot(R_1(:), R_2(:), '.'); axis tight; hold on; fplot(@(x) x, xlim, 'r-');
%         subplot(2,1,2); plot(Rf_1(:), Rf_2(:), '.'); axis tight; hold on; fplot(@(x) x, xlim, 'r-');
%         fprintf('[[cc = %.6f. cc_all = %.6f]]\n', cc_osp, cc_ospf);
%         drawnow;
%         pause(1);
        
        %%
%         [R_1a, Rf_1a, stats1, corticalGain1] = getResponseCorrectedForCorticalState(Gid,cellId, [], 'all');
        
        
%         [R_2, Rf_2, stats2] = getResponseCorrectedForCorticalState(Gid,cellId, [], 'all', Si.corticalGainParams);
        
%         [R_2, Rf_2] = getResponseCorrectedForCorticalState(Gid,cellId, [], 'all', Si.corticalGainParams);
        
        %%
        Si.Gid = Gid;
        Si.cellId = cellId;
        allR_corr_S{i} = Si;
%         Si = getResponseCorrectedForCorticalState(Gid, cellId);
%         allR_corr_S{i} = aa;
%         fprintf('*')
%     catch Merr
%         fprintf('-')
% %         display(Merr);
%     end
%     [L_bin_stimw, R_bin_stimw, windowProfile_best] = getLRbins_windowprofile('stimw', PSTH);
%     S = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_stimw, R_bin_stimw, [], {'osp_gainCorrected'}, redo);
    %%
%     S = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_ph'});
%     S = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, [], {'osp_ph'});
%     if strcmp(class(S), 'MException')
%         S = getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_gainCorrected'}, 1);
%     end
%     getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_gainCorrected_oe'});
%     getOspDataForPsthWindow(Gid, cellId, [], [], L_bin_best, R_bin_best, windowProfile_best, {'osp_gainCorrected_hoe'});
%     progressBar(i);
end


