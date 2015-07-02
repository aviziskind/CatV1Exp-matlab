function [Wcc_pairIdxs, Wcm_pairIdxs, Bcc_pairIdxs, Bcm_pairIdxs] = ...
...% pairsWithinSite_CellCell_S, pairsWithinSite_CellMU_S, pairsBetweenSite_CellCell_S, pairsBetweenSite_CellMU_S
    getCellGroupPairs(cellGroupIds, cellIds, okPairs)    

    nUnits = length(cellIds);
    % cellGroupIds is a vector [1..nCells], containing the group ids of each cell.
    % cellIds is a vector [1..nCells], containing the ids of each cell.
	% cellF1oDCs, nOris, nSps, and nPhs are [1..nGroups] vectors, with the respective parameters of each group

    % 1. Check which groups have the same nori/nsp/nph
%     [uniqueGroupIds, tmp, groupInds] = unique(cellGroupIds);
%     okPairs = false( length(uniqueGroupIds) );    

%     sps_deg = cellfun(@times, sps_pix, num2cell(degPerPixs), 'un', false);
%     sps_deg = sps_pix;

    
%     th = .1;

%     function tf = cmpSpf(sps1, sps2)
%         tf = length(sps1) == length(sps2);
%         return;
%         if isequal(sps1, sps2)
%             tf = true;
%             return;
%         end
%         if length(sps1) ~= length(sps2)
%             tf = false;
%             return;
%         end
%         log1 = log(sps1);
%         log2 = log(sps2);
%         ds = abs((log1 ./ log2) - 1);        
%         tf = all( ds < th);
%     end
    
%     
%     for i = 1:length(uniqueGroupIds)
%        for j = i+1:length(uniqueGroupIds)
%             okPairs(i,j) = cmp1(oris{i}, oris{j}) && cmp1(sps_deg{i}, sps_deg{j}) && cmp1(phs{i}, phs{j});
% %             [nOris(i) nSps(i) nPhs(i)] == [nOris(j) nSps(j) nPhs(j)]);
%             % let bottom-left half of matrix be zero.
% %             okPairs(j,i) = okPairs(i,j);
%        end
%     end
    
    % 2. Generate lists of pairs    
    [pairsWithinCellGroups, pairsBetweenCellGroups] = getPairsWithinAndBetweenGroups(cellGroupIds, okPairs);

    % 3. Split each list by pairs that include a cell multi-cluster and ones that don't   
    W_MU1 = cellIds(pairsWithinCellGroups(:,1)) == 0;   W_MU2 = cellIds(pairsWithinCellGroups(:,2)) == 0;
    B_MU1 = cellIds(pairsBetweenCellGroups(:,1)) == 0;  B_MU2 = cellIds(pairsBetweenCellGroups(:,2)) == 0;
    
    withinGroups_WithMU     =  bitxor( W_MU1, W_MU2 );
    withinGroups_WithoutMU  =  ~ ( W_MU1 | W_MU2 );
    betweenGroups_WithMU    = bitxor( B_MU1, B_MU2 );            
    betweenGroups_WithoutMU =  ~ ( B_MU1 | B_MU2 );            
    
    pairsWcc_V = pairsWithinCellGroups(  withinGroups_WithoutMU, :);
    pairsWcm_V = pairsWithinCellGroups(  withinGroups_WithMU, :);
    pairsBcc_V = pairsBetweenCellGroups( betweenGroups_WithoutMU,:);
    pairsBcm_V = pairsBetweenCellGroups( betweenGroups_WithMU,:);

    mtxSub2ind = @(v) (v(:,2)-1)*nUnits + v(:,1);    
    
    Wcc_pairIdxs = mtxSub2ind(pairsWcc_V);
    Wcm_pairIdxs = mtxSub2ind(pairsWcm_V);
    Bcc_pairIdxs = mtxSub2ind(pairsBcc_V);
    Bcm_pairIdxs = mtxSub2ind(pairsBcm_V);
    

%{    
    Wcc_ori_sp_maxMUs = ori_sp_maxMUs( groupInds( pairsWcc_V(:,1) ) );
    Wcm_ori_sp_maxMUs = ori_sp_maxMUs( groupInds( pairsWcm_V(:,1) ) );
    Bcc_ori_sp_maxMUs = ori_sp_maxMUs( groupInds( pairsBcc_V(:,1) ) ); % even though have different multi-units - just pick 1 of them at random
    Bcm_ori_sp_maxMUs = ori_sp_maxMUs( groupInds( pairsBcm_V(:,1) ) );
    

    Wcc_pairF1oDCs  = zeros(size(pairsWcc_V));
%     Wcc_contrib = zeros(size(pairsWcc_V));
    for i = 1:size(pairsWcc_V,1);
        inds = pairsWcc_V(i,[1:2]);
        Wcc_pairF1oDCs(i,1:2) = cellF1oDCs(inds);
%         n = nCellsInEachGroup(cellGroupIds(i));        
%         Wcc_contrib(i) = iff(any(n == [0, 1]), 0, 2/n);
    end
    
    Wcm_pairF1oDCs    = zeros(size(pairsWcm_V));
    for i = 1:size(pairsWcm_V,1);
        inds = pairsWcm_V(i,[1:2]);
        Wcm_pairF1oDCs(i,1:2) = cellF1oDCs(inds);        
    end
    
    Bcc_pairF1oDCs = zeros(size(pairsBcc_V));
    for i = 1:size(pairsBcc_V,1);
        inds = pairsBcc_V(i,[1:2]);
        Bcc_pairF1oDCs(i,1:2) = cellF1oDCs(inds);
    end
    
    Bcm_pairF1oDCs   = zeros(size(pairsBcm_V));
    for i = 1:size(pairsBcm_V,1);
        inds = pairsBcm_V(i,[1:2]);
        Bcm_pairF1oDCs(i,1:2) = cellF1oDCs(inds);
    end
    
    pairsWithinSite_CellCell_S = struct(...
        'cell1',   num2cell(pairsWcc_V(:,1)),      'cell2',   num2cell(pairsWcc_V(:,2)), ...
        'F1oDC1',  num2cell(Wcc_pairF1oDCs(:,1)),  'F1oDC2',  num2cell(Wcc_pairF1oDCs(:,2)), ...
...%         'contrib', num2cell(Wcc_contribs),  ...
        'ori_sp_max_MU', Wcc_ori_sp_maxMUs);
    pairsWithinSite_CellMU_S   = struct(...
        'cell1', num2cell(pairsWcm_V(:,1)),        'cell2', num2cell(pairsWcm_V(:,2)), ...
        'F1oDC1', num2cell(Wcm_pairF1oDCs(:,1)),   'F1oDC2', num2cell(Wcm_pairF1oDCs(:,2)), ...
...%         'contrib', num2cell(1),...
        'ori_sp_max_MU', Wcm_ori_sp_maxMUs );
    pairsBetweenSite_CellCell_S = struct(...
        'cell1', num2cell(pairsBcc_V(:,1)),        'cell2', num2cell(pairsBcc_V(:,2)), ...
        'F1oDC1', num2cell(Bcc_pairF1oDCs(:,1)),   'F1oDC2', num2cell(Bcc_pairF1oDCs(:,2)), ...
...%         'contrib', num2cell(1),...    
        'ori_sp_max_MU', Bcc_ori_sp_maxMUs );
    pairsBetweenSite_CellMU_S  = struct(...
        'cell1', num2cell(pairsBcm_V(:,1)),         'cell2', num2cell(pairsBcm_V(:,2)), ...
        'F1oDC1', num2cell(Bcm_pairF1oDCs(:,1)),    'F1oDC2', num2cell(Bcm_pairF1oDCs(:,2)), ...
...%         'contrib', num2cell(1),...    
        'ori_sp_max_MU', Bcm_ori_sp_maxMUs );
  %}                        
    fprintf('From %d cells, found : \n   %d within-site cell-cell pairs, %d within-site cell-MU pairs, %d between-site cell cell pairs, and %d between-site cell-MU pairs\n', ...
        length(cellGroupIds), size(pairsWcc_V,1), size(pairsWcm_V,1), size(pairsBcc_V,1), size(pairsBcm_V,1));
    
end


