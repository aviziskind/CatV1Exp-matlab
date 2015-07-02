function displayCellOSP_PSTHs(Gids, glob_idx, allSame)
    nGrps = length(Gids);
    ignoreMultiUnits = true;
    
    S_mcg = load('cellsGroups_movie_fg');
    mcg = S_mcg.movieGroups_fg;    
    [tmp, idx_gids] = intersect([mcg.Gid], Gids);
    
    cellIds = {mcg(idx_gids).cellIds};
    if ignoreMultiUnits
        cellIds = cellfun(@(x) setdiff(x, 0), cellIds, 'un', 0);
    end
    nCells = cellfun(@length, cellIds);
    
    S_data = load('flashedGratingCells_all');
    allCells = S_data.allCells;
    
    [cell_psths, cell_osps] = deal( arrayfun(@(n) cell(1,n), nCells, 'un', 0) );
    
    cell_idxs = cell(1,nGrps);
    for grp_i = 1:nGrps                
        for cell_i = 1:nCells(grp_i)
            Gid = Gids(grp_i); cellId = cellIds{grp_i}(cell_i);
            id = find([allCells.Gid] == Gid  & [allCells.cellId] == cellId, 1);
            if isempty(id)
                error('not found: Gid = %d, celId = %d', Gid, cellId)
            end
            cell_idxs{grp_i}(cell_i) = id;
            cell_psths{grp_i}{cell_i} = allCells(id).PSTH;
            cell_osps{grp_i}{cell_i} = allCells(id).R;
        end
    end
            
    maxNCells = max(nCells);
    
    [idx_matches] = matchUpCellsFromDifferentGroups(cell_osps);
        
    
    r = 5;
    n_row = maxNCells*2*r + (maxNCells-1);
    figure(132); clf;
    for grp_i = 1:nGrps
        for cell_i = 1:nCells(grp_i);
            % psth plot
            cell_i_pos = idx_matches{grp_i}(cell_i);
            
            psth_idx = n_row*(grp_i-1) + (2*r+1)*(cell_i_pos-1)+ [1 : r];
            osp_idx =  psth_idx + r;
            
            h = subplot(nGrps, n_row, psth_idx);
            set(h, 'active', 'outer');
            plotThisPSTH(cell_psths{grp_i}{cell_i})
            if cell_i == 1;
                if (grp_i == 1) && nargin >= 2                    
                    str0 = ['[ ' num2str(glob_idx) ' ]'];
                else
                    str0 = ' ';
                end
                col = 'k';
                if nargin == 3
                    col = iff(allSame, 'k', 'r');
                end                
                st = allCells(cell_idxs{grp_i}(cell_i)).stimType;
                str1 = sprintf('Group # %d', Gids(grp_i)) ;
                idx1 = strfind(st, ':36')+1; idx2 = strfind(st, '('); 
                str2 = st(idx1:idx2-1); str3 = st(idx2:end);                
                allStrs = {str0, str1, str2, str3};
%                 set(gca,'units', 'pixels'); 
                p = get(h, 'position'); 
                scl = @(x, c) x(1)+x(2)*c;
%                 WH = switchh(length(nGrps), {1, 2, 3}, {[.09, .3], [.
                annotation('textbox', [.01, scl(p([2,4]), 0), .09, .2], 'string', allStrs, 'hor', 'center', 'interp', 'none', 'vert','middle', 'color', col, 'fontsize', 9);
%                 text(diff(xlim)*(-1.3), diff(ylim)*.5, allStrs, 'horiz', 'center', 'interpreter', 'none', 'color', col);
%                 set(gca,'units', 'normalized'); 
                
            end
            title(sprintf('Cell # %d', cellIds{grp_i}(cell_i)));

            subplot(nGrps, n_row, osp_idx);
            imageOSP(cell_osps{grp_i}{cell_i}, 'mean:ph', 'OSP', 'nolabels', 'noticks');                        
            3;
            cc_p = allCells(cell_idxs{grp_i}(cell_i)).stats.allWindowStats.cc_p;
            title(sprintf('%.2f', cc_p))
        end
    end
            



end


