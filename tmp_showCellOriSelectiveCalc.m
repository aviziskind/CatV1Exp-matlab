function tmp_showCellOriSelectiveCalc(S)
    
    Gid = S.Gid;
    cellId = S.cellId;
    PSTHdata = S.PSTH;
    
    if isempty(PSTHdata)
        LR_bins = [1 1];
        windowProfile = [];
    else                
        LR_bins = v.PSTH.timeWindow_bins;
        windowProfile = v.PSTH.windowProfile;
    end
    
    redo_flag = 1;  
    getOspDataForPsthWindow(Gid, cellId, [], [], LR_bins(1), LR_bins(2), windowProfile, {'tuningStats'}, redo_flag);
    
    3;



end


%%
S = load('driftingGratingCells_GLFcuw8_degree_all.mat');
idx_ori_cells = arrayfun(@(s) s.cellId > 0 & length(s.ori) > 20, S.allCells);
oriCells = S.allCells(idx_ori_cells);

oriStats = nestedFields(oriCells, 'stats', 'tuningStats', 'oriStats_si');
idx_rep = [oriStats.ori_rep_pval] < .01; % & [oriStats.rsqr] > 0.4;
oriSelStats = nestedFields(oriCells(idx_rep), 'stats', 'tuningStats', 'oriStats_si', 'ori_sel_pval_stats');
p_res = -log10([oriSelStats.resultant_prob]);
p_ttest = -log10([oriSelStats.proj_ttest]);
p_Uall = -log10([oriSelStats.pval_U_allVals]);
p_Tall = -log10([oriSelStats.pval_T_allVals]);

% plot(p_Uall, p_Tall, '.'); drawHorizontalLine(2); drawVerticalLine(2); axis([0 10 0 10]);
% plot(p_Uall, p_res, '.'); drawHorizontalLine(2); drawVerticalLine(2);
% plot(p_Uall, p_ttest, '.'); drawHorizontalLine(2); drawVerticalLine(2);
% plot(p_Uall, p_ttest, '.'); drawHorizontalLine(2); drawVerticalLine(2); axis([0 12 0 12])
% plot(p_Uall, p_res, '.'); drawHorizontalLine(2); drawVerticalLine(2); axis([0 12 0 12])
% plot(p_Uall, p_res, '.'); drawHorizontalLine(2); drawVerticalLine(2);
idx_discrep_U_res = find(xor((p_Uall > 2), (p_res > 2)));
% plot(p_Uall(idx_discrep_U_res), p_res(idx_discrep_U_res), 'ro');
[~, idx_big] = sort( abs( p_Uall(idx_discrep_U_res)-p_res(idx_discrep_U_res) ), 'descend');
oriRepCells = oriCells(idx_rep);
% oriRepCells(473)

% tmp_showCellOriSelectiveCalc( oriRepCells(idx_discrep_U_res(idx_big(1)) ) )
