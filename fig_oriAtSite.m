function fig_oriAtSite

    %%
    S = load('driftingGratingCells_GLFcuw8_degree_SS_ori.mat');
    %%
    allGids = [S.allCells.Gid];
    allCellIds = [S.allCells.cellId];
    
    [Gid_grps, grp_idxs, nUnits_in_group] = uniqueList(allGids);
    CellIds_grp = cellfun(@(idx) allCellIds(idx), grp_idxs, 'un', 0);
    
    hasMU = cellfun(@(cids) any(cids == 0), CellIds_grp);
    nCellsInGroup = cellfun(@(cids) nnz(cids > 0), CellIds_grp);
    
    minNCells = 3;
    
    oriStats = nestedFields(S.allCells, 'tuningStats', 'oriStats_ss', 1);
    dMU = [oriStats.Dori_pref_smlSpkMU];
    has1Outlier = cellfun(@(idxs) nnz(dMU(idxs) > 45 == 1), grp_idxs);
    
    idx_grp_ok = find(hasMU & (nCellsInGroup >= minNCells) & has1Outlier);
    nOK = nnz(idx_grp_ok);
    fig_id = 110;
    
    %%
%     for grp_i = 1:length(idx_grp_ok)
        grp_i = 11;
        oriMax = 360;
        smoothW = 1;
        figure(fig_id); clf; hold on; box on;
        
        grp_idx = idx_grp_ok(grp_i);
        Gid = Gid_grps(grp_idx);  % use Gid = 5144
        cellIds_orig = CellIds_grp{grp_idx};
        cell_idxs = grp_idxs{grp_idx};
       
        
        idx_cell_use = 1:length(cell_idxs);
        idx_cell_use(2) = [];
        cell_idxs = cell_idxs(idx_cell_use);
        
        cellIds = cellIds_orig(idx_cell_use);
        
        oris = [S.allCells(cell_idxs(1)).ori(:)];
        locData = [S.allCells(cell_idxs(1)).locData];
        
        dx = diff(oris(1:2));
        ext = @(x) [x; x(end)+dx];
        wrp = @(x) [x; x(1)];

        cols = ['kbrcgmy'];
        h = [];
        smooth_w = .6;
        shift = 0;
        nShift = round(shift / dx);
        oris_fine = linspace(oris(1), oris(end) + diff(oris(1:2)), length(oris)*10+1);
        nCells = nnz(cellIds);
        allBarLayers = {};
        
        for ci = 1:length(cellIds)
%             cellId_show = find(cellIds(ci)==cellIds,1);
            stat_i = oriStats(cell_idxs(ci));
            dir_pref = stat_i.dir_pref_deg + shift;
            line_s = iff(cellIds(ci) == 0, ':', ':');
            col = iff(cellIds(ci) == 0, 'k', cols(ci));
%             oris_shft = mod(oris, 360);
            
            r_k_dir = stat_i.r_k_dir;
            if smooth_w > 0
                r_k_dir = gaussSmooth(r_k_dir, smooth_w, [], 1);
            end
            
            otc = circshift(r_k_dir, nShift);
            h(ci) = plot(ext(oris), wrp(otc), [col, 'o' line_s], 'markersize', 4, 'linewidth', 2, 'linestyle', ':');
            
            p = stat_i.oriParams;
            if cellIds(ci) > 0
                otc_i = gaussOri(p.A, p.sigma, p.B, circDist(oris_fine,dir_pref));
                
                oris_fine_i = dir_pref + linspace(-90, 90, length(oris)*10+1);
                otc_i = gaussOri(p.A, p.sigma, p.B, oris_fine_i-dir_pref);

                [oris_fine_i_wrp, otc_i_wrp] = wrapCurve(oris_fine_i, otc_i, 360);
                plot(oris_fine_i_wrp, otc_i_wrp, [col '-'], 'linewidth', 2);
            end
            
            plotWidths = 1;
            showOriLocal = 1;
            showOriGlobal = 0;
            showPrefOri = 1;
            showJackKnifeStdErr = 0;
            if plotWidths && cellIds(ci) > 0
                ori_w_local = stat_i.w_ori_local;
                ori_w_global = stat_i.w_ori_global;
                
%                 fprintf('%.1f, %.1f\n', ori_w_local, ori_w_global);
                fprintf('Cell %d (id = %d), pref = %.1f. local w = %.1f, global w = %.1f\n', ci, cellIds(ci), dir_pref, ori_w_local, ori_w_global);
%                 width_bars_y_range = [20, 22];
                width_bars_y_range = [65, 85];
                bar_y_spacing = 3;
                
                ori_w_local_wasnan = isnan(ori_w_local);
                if ori_w_local_wasnan
                    ori_w_local = stat_i.orig.w_ori_local;
                end
                ext_factor = 1.1;
                x_pos = dir_pref;
                barInterval = sort( x_pos + ori_w_local*[-1, 1]*ext_factor );
                    
                curBarLayer = 1;
                foundSpace = false;
                while ~foundSpace
                    
                    anyOvlp = false;
                    if length(allBarLayers) < curBarLayer
                        allBarLayers{curBarLayer} = {};
                    end
                    if ~isempty(allBarLayers{curBarLayer})
                        for i = 1:length(allBarLayers{curBarLayer})
                            if intervalsOverlap(allBarLayers{curBarLayer}{i}, barInterval)
                                anyOvlp = true;
                            end
                        end
                    end
                    
                    if ~anyOvlp
                        foundSpace = true;
                    else
                        curBarLayer = curBarLayer + 1;
                    end
                    
                end 
                allBarLayers{curBarLayer}{end+1} = [barInterval];
                
                
%                 y_pos = diff(width_bars_y_range)*(ci-1)/nCells + width_bars_y_range(1);
                y_pos = width_bars_y_range(1) + (curBarLayer-1) * bar_y_spacing;
                if showPrefOri
                    
                    prefHeight = 1.5;
                    plot(x_pos*[1,1], y_pos+[-1,1]*prefHeight, [col, '-'], 'linewidth', 2);
                    
                end
                
                
                
                if showOriLocal
                    %%
                    
                   
                    linestyle = iff(ori_w_local_wasnan, ':', '-');                        
                    plot(x_pos-ori_w_local*[-1, 1], y_pos*[1,1], col, 'linewidth', 3, 'linestyle', linestyle);
                    
                 
                    
                end
                   if showJackKnifeStdErr
                        plot(x_pos-stat_i.error_jack.ori_pref *[-1, 1], y_pos*[1,1], 'k', 'linewidth', 2);
                   end
                    
                if showOriGlobal
                    plot(x_pos-ori_w_global*[-1, 1], y_pos*[1,1], [col, '-'], 'linewidth', 2);
                    endHeight = 1.5;
                    plot(x_pos-ori_w_global*[1,1], y_pos+[-1,1]*endHeight, [col, '-'], 'linewidth', 2);
                    plot(x_pos+ori_w_global*[1,1], y_pos+[-1,1]*endHeight, [col, '-'], 'linewidth', 2);
                end
              
                
            end
            
            
        end
        
        %%
        label_fontSize = 13;
        legend_fontsize = 10;
        axis_fontSize = 11;
        
        for ci = 1:length(cellIds)
            stat_i = oriStats(cell_idxs(ci));
            dir_pref = stat_i.dir_pref_deg + shift;
            col = iff(cellIds(ci) == 0, 'k', color_s(ci-1));
            if cellIds(ci) == 0
               drawVerticalLine(dir_pref, 'color', cols(ci), 'linewidth', 2, 'linestyle', ':') 
            end
        end
        set(gca, 'xtick', [0:45:360], 'xlim', [0 360])
        legend(h, ['Multiunit'; legendarray('Cell ', 1:nnz(cellIds))], 'location', 'best', 'fontsize', legend_fontsize)
        xlabel('Direction (degrees)', 'fontsize', label_fontSize);
        ylabel('Mean Firing Rate (Hz)', 'fontsize', label_fontSize);
%         title(sprintf('Orientation/Direction tuning curves (site # %d)', locData.LocId))
        
        
%     end
    
    %%
        figureFolder = [CatV1Path 'Figures' filesep 'DegreePaper' filesep];
        fig_filename = sprintf('%sFigure1_orisAtASite.pdf', figureFolder);
        set(fig_id, 'color', 'w', 'windowStyle', 'normal', 'position', [1100   400   780   380]);
        set(gca, 'fontSize', axis_fontSize)
        %%
        export_fig(fig_id, 'cymk', fig_filename);
        3;
%%
%         export_fig(60, [CatV1Path 'test\test3.pdf']);
    
%%

end


        
        
function tf = intervalsOverlap(int1, int2)
    tf = ~ ((int2(2) < int1(1)) || (int2(1) > int1(2)) );

end