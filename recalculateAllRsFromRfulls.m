% function recalculateAllRsFromRfulls

    [gratingTypeId, gratingType] = curGratingType;
    fprintf(['Recalculating R (from R_full) for ' gratingType ' gratings ...\n'])
    
    s = who('celldata__*');
    [Gids] = cellfun(@getGidCellIdFromVarname, s);
    cell_gTypeId = flashedOrDrifting(Gids);    
    
    idx = find(cell_gTypeId == gratingTypeId);
        
    nCells = length(idx);
    fprintf('Attempting to recalculate R for %d cells ... \n', nCells);
    nRecalced = 0;
    progressBar('init-', nCells, 60);
    for i = 1:nCells
        progressBar;
        varname = s{idx(i)};
        S = eval(varname);
        if isfield(S, 'OSP') && isfield(S.OSP, 'stats') && isfield(S.OSP, 'R_full') && ~iscell(S.OSP.R)
            Gid = S.Gid;            
            assert(flashedOrDrifting(Gid) == gratingTypeId)
            [uOri, uSp, uPh, nPermsOrCycles, nReps] = dbGetUniqueOriSpPh('Gid', Gid);
            nRepsTotal = nPermsOrCycles * nReps;
            OSP = S.OSP;
            R_old = OSP.R;
            
            R_full = S.OSP.R_full;
            f = 1;
            if isfield(OSP, 'R_full_factor')
                f = OSP.R_full_factor;
            end
            if ~isstruct(f), f = struct('scale', f, 'idxnan', []); end
            
            if strcmp(class(R_full), 'uint8')
                R_full = single(R_full)/f.scale;
                R_full(f.idxnan) = nan;                
            end            
            [nOri, nSp, nPh, nTrialsMax] = size(R_full);
            nTrialsTheseParams = sum(~isnan(R_full),4);            
            if (flashedOrDrifting(Gid) == 2)  % for drifting gratings - have to discard first cycle which may contain artifacts (of sudden change from no stimulus)
                nCyclesToSkip = 1;
                nCycles = ceil(nPermsOrCycles);
                while (nCyclesToSkip >= nCycles-1)
                    nCyclesToSkip = nCyclesToSkip - 1;
                end
                cycleIdx = 1+nCyclesToSkip:nCycles;
%                 cycleIdx = nCycles;
                okTrialsInds = bsxfun(@plus, [0:nReps-1]'*nCycles, cycleIdx)';
                okTrialsInds = okTrialsInds(:)';
            else
                okTrialsInds = 1:nTrialsMax;
            end
            R = nanmean(R_full(:,:,:,okTrialsInds), 4);
            R(nTrialsTheseParams == 0) = 0;
            
            assert(~any(isnan(R(:))));
            
            a = isequal(R, R_old);
            
            OSP.R = R;        
            nRecalced = nRecalced + 1;
            eval([varname '.OSP = OSP;';]);       
        end
    end
    progressBar('done');
    fprintf('Recalculated stats for %d cells ... \n', nRecalced);

    


