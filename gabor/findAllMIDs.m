function findAllMIDs(Gids, cellIds)

    global CPU_id
    redo = 0;

    midDir = [experimentDBPath 'MID' filesep];
    if exist([midDir 'stop_running'], 'file')
        movefile([midDir 'stop_running'], [midDir 'stop_runnin'])
    end

     if nargin < 2
        do128 = 1;
        S = load('usedCells.mat');
        if do128
            idx_start = 1;
        else
            idx_start = find(S.usedGids == 1822, 1, 'last')+1;
        end        
        idx = idx_start:length(S.usedCellIds);
        Gids = S.usedGids(idx);
        cellIds = S.usedCellIds(idx);
        
%         [Gids, cellIds] = getAllGids('f');        
    end
    nCells = length(cellIds);

%     doMID_reps = {1, 2, 'all'}; 
    doMID_reps = {1}; 
    downSampFactors = [1];
    trialModes = {'all', 'odd', 'even'};
%     trialModes = {'odd', 'even'};
%     timeWindows = {[29, 62], [58, 92], 'best'};
%     timeWindows = {'best'};
    timeWindows = {[29, 62], [58, 91]};
    responseType = curResponseType('');
    
    haveMultipleSessions = ~isempty(CPU_id);
    CPU_str = iff(isempty(CPU_id), '', sprintf('(CPU_%d) ', CPU_id));
    CPU_s = iff(isempty(CPU_id), '', sprintf('%d', CPU_id));
    
    nnreps = length(doMID_reps);
    nnsamp = length(downSampFactors);
    nWindows = length(timeWindows);
    nTrialModes = length(trialModes);
    
%     downSampFactor = 2;
    ndone = 0;
%     setMatlabWindowTitle(sprintf('%sscanning...', CPU_str));
    for cell_i = 1:nCells
        Gid = Gids(cell_i);
        cellId = cellIds(cell_i);

        if exist([midDir 'stop_running'], 'file')
            return;
        end        

        allFileNames = cell(nnreps, nnsamp, nWindows, nTrialModes);
        for di = 1:nnsamp
            for ri = 1:nnreps
                nrep = doMID_reps{ri};
                if isnumeric(nrep), 
                    nrep = sprintf('%drep', nrep);
                end
                
                for wi = 1:nWindows
                    timeWindow = timeWindows{wi};
                    for ti = 1:nTrialModes
                        trialMode = trialModes{ti};
                        if (strcmp(nrep, '1rep')) || strcmp(trialMode, 'all');
                            allFileNames{ri, di, wi, ti} = getName('MID_file', Gid, cellId, downSampFactors(di), nrep, timeWindow, trialMode, responseType);
                        else
                            allFileNames{ri, di, wi, ti} = nan;
                        end
                    end
                end
            end
        end            
        
        filesPresent = cellfun(@(s) any(isnan(s)) || exist(s, 'file'), allFileNames);
                
        anythingToDo = any(~filesPresent(:)) || redo;
                   
        if ~anythingToDo
            fprintf('[%d:%d] ', Gid, cellId);
            ndone = ndone+1;
            if ndone > 11
                ndone = 0;
                fprintf('\n');
            end
            continue;
        end

        tryGetLock = nCells > 1;
        
        if tryGetLock
            lock_name = sprintf('MID_%d_of_%d__Group_%d_cell_%d.tmp', cell_i, nCells, Gid, cellId);
%             lock_name = sprintf('%s_workingOn_%d_of_%d__Group_%d_cell_%d.tmp', midDir, cell_i, nCells, Gid, cellId);
            gotLock = lock_createLock(lock_name);
            
        end
        
        if tryGetLock && ~gotLock
            continue;
        end
                    
            
        windowTitle = sprintf('%sDoing Grp %d, cell %d', CPU_str, Gid, cellId);
        setMatlabWindowTitle(windowTitle);
        
        fprintf('*********** Processing Group %d Cell %d (%s) ********** \n', Gid, cellId, outOf(cell_i, nCells));
        
        for di = 1:nnsamp
            for ri = 1:nnreps
                nrep = doMID_reps{ri};
                if isnumeric(nrep), 
                    nrep = sprintf('%drep', nrep);
                end
                
                for wi = 1:nWindows
                    timeWindow = timeWindows{wi};
                    for ti = 1:nTrialModes
                        trialMode = trialModes{ti};
                        if ~filesPresent(di, ri, wi, ti)
%                         if (strcmp(nrep, '1rep')) || strcmp(trialMode, 'all');
                            opt = struct('frameMode', nrep, 'downSampFactor',downSampFactors(di), 'timeWindow', timeWindow, 'trialMode', trialMode, 'responseType', responseType);
                            mid_findMostInfoDim(Gid, cellId, opt);                                
%                         end
                        end
                    end
                end
            end
        end            
        
                                       
%         catch MErr
%             if stopIfError
%                 rethrow(MErr);
%             else
%                 fprintf('Error encountered while doing Gid = %d\n', Gid);           
%                 fprintf('* Identifier: %s\n* Message: %s\n', MErr.identifier, MErr.message);        
%                 errorGids = [errorGids, Gid]; %#ok<AGROW>
%             end
%             
%         end
        if tryGetLock
            lock_removeLock(lock_name);
        end

%         if (nCells > 1) && exist(tmpFileName_full, 'file');
%             delete(tmpFileName_full);
%             % if are working in parallel, and this is the last datafile to
%             % be completed, and we're on the workstation, send a message that we're done.                        
%             if strcmp(getenv('computername'), 'AVI-WORK-PC') && (cell_i == nCells)
%                 s = dir([midDir '_workingOnGroup*']);
%                 if isempty(s)
%                     disp('Sending Email to self...');
% %                     sendEmailToSelf('Completed last datafile in the list');
%                 end
%             end
%             
%         end
    end
   
    setMatlabWindowTitle('Done');
%     if ~isempty(errorGids)
%        fprintf('Encountered errors with the following Gids : ');
%        disp(errorGids);
%     end

       
end