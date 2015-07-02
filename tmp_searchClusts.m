function tmp_searchClusts    
    % search for the Group in Supplementary Figure 3
    
    S = load('cellsGroups_grating_dSf_DB');
    Gids = [S.gratingGroups_dSf.Gid];
    NC = 4;
    idx_n = find( arrayfun(@(s) nnz( s.cellIds > 0)==NC, S.gratingGroups_dSf) );
    Gids_n = Gids(idx_n);
    
%     S = load('cellsGroups_GLFcuw8_grating_dOr');        
    
    S_indiv = load('indivCells_DB_grating_dSf');
%     idx0 = 19+43*[0,1,2,3];
    
    nCells_ok_n = zeros(1, length(idx_n));
    F1oDCs = zeros(4, length(idx_n));
    
    for i = 1:length(idx_n)
        grp_data = S.gratingGroups_dSf(idx_n(i));
        Gid = grp_data.Gid;        
        cellIds = grp_data.cellIds;
        cellIds = cellIds(cellIds > 0);
        
        nSpf = length(grp_data.spPeriod_pix);           
        grpF1oDCs = zeros(NC,1);
        
        for j = 1:length(cellIds)
            vn = getName('celldata', Gid, cellIds(j));
            s = S_indiv.(vn);
            R = s.OSP.R;
%             [~, idx_best] = maxElement( mean(R,3) );
%             cell_ok(j) = s.OSP.stats.tuningStats.spfStats_si.cellOK;            
            
            phs = s.OSP.ph;
            [nOri, nSp, nPh] = size(R); %#ok<NASGU>
            meanR = mean(R, 3);                                                

            [maxR_av, oriSp_maxR_avP] = maxElement(meanR);
            [ori_i_av, sp_i_av]      = elements(oriSp_maxR_avP);                
            phaseTC_atMax_avP = squeeze(R(ori_i_av, sp_i_av, :));                    
            F1oDC_maxR_avP = getF1oDC(phs, phaseTC_atMax_avP, 360 ); % phase_max is always 360 for phases.
            
            grpF1oDCs(j) = F1oDC_maxR_avP;
        end
        F1oDCs(:,i) = grpF1oDCs;
        
        nCells_ok_n(i) = nnz(grpF1oDCs > 1) == 2;
        
        if nCells_ok_n(i)
            fprintf('Gid = %d : %.2f, %.2f, %.2f, %.2f\n', Gid, grpF1oDCs );
%             spfs_cpd = 1./(grp_data.spPeriod_pix * grp_data.degPerPix)
%             figure(5); clf; 
%             for k = 1:5
%                 subplot(5,1,k);
%                 plot(spfs_cpd, stcs(k,:), '.-')
%             end
%             
%             xlabel(sprintf('Gid = %d', Gid))
            3;
        end
        
    end
    G_ncells = Gids_n(find(nCells_ok_n))
    
    3;
end




%{
function tmp_searchClusts    
    % search for the Group in Supplementary Figure 1
    
    S = load('cellsGroups_grating_dOr_DB');
    Gids = [S.gratingGroups_dOr.Gid];
    idx5 = find( arrayfun(@(s) nnz( s.cellIds > 0)>=5, S.gratingGroups_dOr) );
    Gids5 = Gids(idx5);
    
%     S = load('cellsGroups_GLFcuw8_grating_dOr');        
    
    S_indiv = load('indivCells_DB_grating_dOr');
%     idx0 = 19+43*[0,1,2,3];
    
    nCells_ok3 = zeros(1, length(idx5));

    for i = 1:length(Gids5)
        Gid = Gids5(i);
        group_idx = find([S.gratingGroups_dOr.Gid] == Gid, 1);
        grp_data = S.gratingGroups_dOr(group_idx);
        assert(grp_data.Gid == Gid);
        cellIds = grp_data.cellIds;
        cellIds = cellIds(cellIds > 0);
        cell_ok = zeros(1, length(cellIds));
        nOri = length(grp_data.ori_deg);           
        otcs = zeros(5, nOri); 
        for j = 1:length(cellIds)
            vn = getName('celldata', Gid, cellIds(j));
            s = S_indiv.(vn);
            cell_ok(j) = s.OSP.stats.tuningStats.oriStats_si.cellOK;            
            otcs(j,:) = mean(s.OSP.R, 3)';            
        end
        nCells_ok3(i) = nnz(cell_ok) == 3;
        
        if nnz(cell_ok) == 3
            figure(5); clf; 
            for k = 1:5
                subplot(5,1,k);
                plot(otcs(k,:))
            end
            
            xlabel(sprintf('Gid = %d', Gid))
            3;
        end
        
    end
    G_ncells = Gids5(find(nCells_ok3))
    
end

%}

%   447
%         1528
%         1715




% function tmp_searchClusts
%     global GidsAfter5236
%     Gids = getAllGids('o');
%     
%     couldBe = zeros(1, length(Gids));
%     for i = 1:length(Gids)
%         clustIds = getClusterCellAssignment(Gids(i));       
%         if isempty(clustIds) || ~iscell(clustIds)
%             continue;
%         end
%         clustIds = clustIds(~cellfun(@isempty, clustIds));
%         clustIds1 = cellfun(@(x) x(1), clustIds);
%         
%         if any(clustIds1 == 2) && any(clustIds1 == 8)
%             couldBe(i) = 1;
%         end        
%     end
%         
%     candidates = Gids( find(couldBe) )
% 
%     3;
% 
% end

