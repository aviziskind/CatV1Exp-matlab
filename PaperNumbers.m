
%% generate
regenerate = 0;
if regenerate 
    curGratingType(1);
    curCmpType('phase'); generateGratingCellsDatafile;
    curCmpType('degree'); generateGratingCellsDatafile;

    curGratingType(2);
    curCmpType('phase'); generateGratingCellsDatafile;
    curCmpType('degree'); generateGratingCellsDatafile;
end

%% Load Cells

    S_degree_do = load('driftingGratingCells_GLFcuw8_degree_ori');
    S_degree_ds = load('driftingGratingCells_GLFcuw8_degree_spf');
    S_degree_fo = load('flashedGratingCells_GLFcuw8_degree_ori');
    S_degree_fs = load('flashedGratingCells_GLFcuw8_degree_spf');
    S_phase_d = load('driftingGratingCells_GLFcuw8_phase.mat');
    S_phase_f = load('flashedGratingCells_GLFcuw8_phase.mat');
    
    S_labels = {'degree-drifting-ori', 'degree-drifting-spf', 'degree-flashed-ori', 'degree-flashed-spf', 'phase-drifting', 'phase-flashed'};
    S_cmp = {'degree', 'degree', 'degree', 'degree', 'phase', 'phase'};
    
findAvFiringRate = 0;
if findAvFiringRate 
    %% find average firing rate of drifting grating cells to best response
    idx_cell = arrayfun(@(s) s.cellId > 0, S_degree_ds.allCells);
    allCells = S_degree_ds.allCells(idx_cell);
    N = length(allCells);
    fr_max = zeros(1,N);
    for i = 1:N
        meanR = mean(allCells(i).R, 3);
        fr_max(i) = max(meanR(:));        
    end
end        
    
compareBnBlW = 0;
if compareBnBlW
    %% Comparison of Bn, Bl, w
    idx_spf = find(arrayfun(@(s) strncmp(s.stimType, 'Grating:Spat', 10), S_degree_ds.allCells) & ([S_degree_ds.allCells.cellId]' > 0) );
    spfCells = S_degree_ds.allCells(idx_spf);
    spfStats = nestedFields(spfCells, 'stats', 'tuningStats', 'spfStats_si');    
    params = [spfStats.SLNparams];
    w = [params.w]';
    w_spf = [spfStats.w_spf]';
    Bn = [spfStats.Bn]';

    c_Bn_wspf = corr(w_spf, Bn);
    p_Bn_wspf = polyfit(w_spf, Bn, 1);
    fprintf('Bn vs w_spf : best fit: Bn = %.3f w_spf + %.3f.  cc = %.3f\n', p_Bn_wspf(1), p_Bn_wspf(2), c_Bn_wspf);

    c_w_wspf = corr(w_spf, w);
    p_w_wspf = polyfit(w_spf, w, 1);
    fprintf('w  vs w_spf : best fit: w = %.3f w_spf + %.3f.  cc = %.3f\n', p_w_wspf(1), p_w_wspf(2), c_w_wspf);
end
    
%% Collect
% all_S = {S_degree_do.allCells, S_degree_ds.allCells, S_degree_fo.allCells, S_degree_fs.allCells, S_phase_d.allCells, S_phase_f.allCells};
allCells_C = {S_degree_do.allCells, S_degree_ds.allCells, S_degree_fo.allCells, S_degree_fs.allCells, S_phase_d.allCells, S_phase_f.allCells};

% allCells_deg_C = {S_degree_do.allCells, S_degree_ds.allCells, S_degree_fo.allCells, S_degree_fs.allCells};
% allCells_ph_C = {S_phase_d.allCells, S_phase_f.allCells};


%% how many cells in total did we use?

% GC_C_ph = cellfun(@(allCl) [allCl.Gid] * 1000 + [allCl.cellId], allCells_ph_C, 'un', 0);
% GC_C_deg = cellfun(@(allCl) [allCl.Gid] * 1000 + [allCl.cellId], allCells_deg_C, 'un', 0);
% 
% GC_ph = unique([GC_C_ph{:}]);
% 



%% Get Location array, add site info
NN = length(allCells_C);
all_loc_C = cell(1,NN);
for i = 1:NN
    all_loc_C{i} = [allCells_C{i}.locData];
    for j = 1:length(allCells_C{i})
        Gid = allCells_C{i}(j).Gid;
        sd = siteDataFor(Gid);
        date_str = sd.dataFileInfo.dateCreated;
        
        id_col = strfind(allCells_C{i}(j).stimType, ':');
        stimType = allCells_C{i}(j).stimType(1:id_col(2)-1);        
        
        stimType_short = switchh(stimType, {'Grating:Orientation Batch', 'Grating:Spatial Frequency Batch',  'Movie:Flashed_Gratings'}, {'orient', 'spat', 'flashed'});
        all_loc_C{i}(j).date = date_str;
        all_loc_C{i}(j).stimType = stimType_short;
        all_loc_C{i}(j).loc_stim = sprintf('%d_%s', all_loc_C{i}(j).LocId, stimType_short);
        all_loc_C{i}(j).Gid = Gid;
        all_loc_C{i}(j).cellId = allCells_C{i}(j).cellId;
        all_loc_C{i}(j).GC = allCells_C{i}(j).Gid * 10000 + allCells_C{i}(j).cellId;
        all_loc_C{i}(j).cmpType = S_cmp{i};
    end    
end
all_loc = [all_loc_C{:}];
%% remove duplicates 
% some cells are repeated (appear in degree & phase tuning)
all_GC = [all_loc.GC];
[uGC, GC_list] = uniqueList(all_GC);

all_cells = repmat(blankStruct(all_loc(1)), 1, length(uGC));
for i = 1:length(uGC)
    all_cells(i) = all_loc(GC_list{i}(1));
    newCmpType = cellstr2csslist( {all_loc(GC_list{i}).cmpType} );
    all_cells(i).cmpType = newCmpType;
end

%% Does every group have a multi-unit? (Yes)
uGids = unique([all_loc.Gid]); 
MU_idx = find([all_loc.cellId] == 0);
MU_Gids = [all_loc(MU_idx).Gid]; %#ok<FNDSB>
Gids_with_no_MU = setdiff(uGids, MU_Gids);   % empty matrix
fprintf('%d groups have no multi-unit\n', length(Gids_with_no_MU));

%% Just multiunits of each site
MU_idx = find([all_cells.cellId] == 0);
all_cell_MU = all_cells(MU_idx);


%% How many orientations were followed by spat freq/flashed?
[uLocs, locList] = uniqueList([all_cell_MU.LocId]);
loc_str = cell(1,length(uLocs));
for i = 1:length(uLocs)    
    Gids_here = [all_cell_MU(locList{i}).Gid];
   idx = ord( {all_cell_MU(locList{i}).date} );    
   loc_str{i} = cellstr2csslist(  sort({all_cell_MU(locList{i}(idx)).stimType}) );    
end
varBreakdown(loc_str)

find(strcmp('spat, orient', loc_str))


%% Display Total # Cats/Pens/Groups/Cells for each set of experiments.
for i = 1:4
    loc = all_loc_C{i};
    idx_cells = [loc.cellId] > 0;
    loc = loc(idx_cells);
    nCats = length(unique([loc.CatId]));
    nPens = length(unique([loc.PenId]));
    nLocs = length(unique([loc.LocId]));
    nGrps = length(unique([loc.Gid]));    
    nCells = length(unique([loc.GC]));    
    
    fprintf('For %s : %d Cats, %d Penetrations, %d Locations, %d Groups, %d Cells, \n', ...
        S_labels{i}, nCats, nPens, nLocs, nGrps, nCells);
    

    
end


%% display total number of cells used, accounting for duplicate experiments at a recording location
% nCells =



%%
% how old each of the cats were (verify 3-6 months).
[date_birth_str, date_exp_str, lab_name] = getFieldsFromDatabaseTable(hnd, {'DTM_BIRTH', 'DTM_EXPERIMENT', 'TXT_LAB_NAME'}, 'TBL_ANIMALS');
idx_use = find(strncmp(lab_name, 'K', 1));
date_birth_str = date_birth_str(idx_use); date_exp_str = date_exp_str(idx_use); lab_name(idx_use)
date_birth = cellfun(@datenum, date_birth_str);
date_exp   = cellfun(@datenum, date_exp_str);
ndays_old = (date_exp-date_birth);
nmonths_old = ndays_old/30;





    %%
    
    %% d
    deg_loc = [deg_d_loc, deg_f_loc];
    deg_d_loc = [S_degree_d.allCells.locData];
    deg_f_loc = [S_degree_f.allCells.locData];
    phase_d_loc = [S_phase_d.allCells.locData];
    phase_f_loc = [S_phase_f.allCells.locData];

    
    f
    
    all_loc = [deg_d_loc, deg_f_loc, phase_d_loc, phase_d_loc];
    a = 3;
    
%     fprintf('Total of 

    
    
    
    
    
