function [gids, cellIds] = getAllGids(gratingTypes, sortBy)
    if (nargin < 1) || isempty(gratingTypes);
        gratingTypes = 'fd';
    end
    if (nargin < 2) || isempty(sortBy);
        sortBy = 'Gid';
    end
    
    gids = [];
    nspk = [];
    cellIds = {};
    doFlashedGratings = any(gratingTypes == 'f');
    doF16 = ~isempty(strfind(gratingTypes, 'f16'));
    doF41 = ~isempty(strfind(gratingTypes, 'f41'));
    doF100 = ~isempty(strfind(gratingTypes, 'f100'));
    doAllFlashedGratings = ~doF16 && ~doF41 && ~doF100;
    
    doAllDriftingGratings = any(gratingTypes == 'd');
    doSpatDriftingGratings = any(gratingTypes == 's');
    doOriDriftingGratings = any(gratingTypes == 'o');
    outputCellIds = nargout > 1;
    
    if doFlashedGratings        
        S_f = load(getFileName('Groups', 'movie_fg'));
        gids = [S_f.movieGroups_fg.Gid];
        nspk = arrayfun(@(grp) sum(grp.nSpikes), S_f.movieGroups_fg);
        if outputCellIds
            cellIds = {S_f.movieGroups_fg.cellIds};            
        end
        if ~doAllFlashedGratings
            frameLength_ms = floor([S_f.movieGroups_fg.frameLength_ms]);
            idx = false(1,length(gids));
            idx(frameLength_ms == 16) = doF16;
            idx(frameLength_ms == 41) = doF41;
            idx(frameLength_ms == 100) = doF100;
                
            gids = gids(idx);
            nspk = nspk(idx);
            if outputCellIds
                cellIds = cellIds(idx);
            end
        end
        
        
    end
    if doAllDriftingGratings || doSpatDriftingGratings
        S_dSf = load(getFileName('Groups', 'grating_dSf'));  
        dg_gids = [S_dSf.gratingGroups_dSf.Gid];
        dg_nspk = arrayfun(@(grp) sum(grp.nSpikes), S_dSf.gratingGroups_dSf);
        gids = [gids(:); dg_gids(:)];
        nspk = [nspk(:); dg_nspk(:)];        
        if outputCellIds
            cellIds = [cellIds, {S_dSf.gratingGroups_dSf.cellIds}];            
        end        
    end
    if doAllDriftingGratings || doOriDriftingGratings
        S_dOr = load(getFileName('Groups', 'grating_dOr'));  
        dg_gids = [S_dOr.gratingGroups_dOr.Gid];
        dg_nspk = arrayfun(@(grp) sum(grp.nSpikes), S_dOr.gratingGroups_dOr);
        gids = [gids(:); dg_gids(:)];
        nspk = [nspk(:); dg_nspk(:)];        
        if outputCellIds
            cellIds = [cellIds, {S_dOr.gratingGroups_dOr.cellIds}];            
        end                
    end
    
    
    switch sortBy
        case 'Gid', % not sorted by Gid if combine of the types above
            idx = ord(gids, 'ascend');
            gids = gids(idx);
            if outputCellIds
                cellIds = cellIds(idx);            
            end
        case 'nSpikes',
            idx = ord(nspk, 'ascend');
            gids = gids(idx);
            if outputCellIds
                cellIds = cellIds(idx);
            end
    end

    if outputCellIds
        cellIds = cellfun(@(x) x(x>0), cellIds, 'un', 0);
        nCellsPerGrp = cellfun(@length, cellIds);
        cellIds = double([cellIds{:}]);
        gids_C = arrayfun(@(G,n) G(ones(1,n)), gids(:), nCellsPerGrp(:), 'un', 0);
        gids = [gids_C{:}];        
    end
    gids = gids(:);
    
%     gids = setdiff(gids, 4522);
end