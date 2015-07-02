function spikes = getSpikes(Gid, cellId, outputType, additionalFeatures, matchSpiker, file_opt)
    % Can return either all spikes for a group of cells, or just 1 cell.
    %       groupSpikes = dbGetSpikes(Gid, [])  % all cells in the group
    %       groupSpikes = dbGetSpikes(Gid, cellId)  % only this cell
    %       groupSpikes = dbGetSpikes(..., outputType) %outputType = 'ms'/'sec'/'tick' 
    %       groupSpikes = dbGetSpikes(..., outputType, keepTetrodeDataFlag)
    
    %%% GET GROUP SPIKE DATA

    if nargin < 2
        cellId = [];
    end
    if nargin < 3
        outputType = [];
    end
    if nargin < 4
        additionalFeatures = [];
    end
    if nargin < 5
        matchSpiker = curMatchDB;
    end
    if nargin < 6
        file_opt = struct;
    end
            
    if matchSpiker
        spikes = dbGetSpikes(Gid, cellId, outputType, additionalFeatures, matchSpiker);
        return;
    end
        
    
    if isscalar(Gid)
        propFileName = getFileName('properties', Gid, [], file_opt);
        S_prop = load(propFileName);
        spikeTimes = double( S_prop.position );
        nSpk = length(spikeTimes);
        % add ability to select  additionalFeatures        
        if ~isempty(additionalFeatures)            
            if ischar(additionalFeatures)
                additionalFeatures = {additionalFeatures};
            elseif (additionalFeatures == 1)
                additionalFeatures = {'negAmps'};
            end
            nFet = length(additionalFeatures);
            fet = cell(1, nFet);
            for fi = 1:nFet
                fet{fi} = double( S_prop.(additionalFeatures{fi}) );
            end
            fet = [fet{:}];
        else
            fet = zeros(nSpk, 0);
        end
                        
        
%         cluFileName = getFileName(curGroupingType(''), Gid);
%         S_clu = load(cluFileName);
        spkCellIds = double( getCellSorting(Gid, 'cells', [], file_opt) )';        
               
        
        assert(length(spkCellIds) == nSpk);
    else
        groupSpikes = Gid;
        spikeTimes = groupSpikes(:,1);
        spkCellIds = groupSpikes(:,2);        
        fet = groupSpikes(:,3:end);
        
    end
        
        
    %%% OUTPUT  GROUP/CELL DATA
        
    % 2. Select which spikes to keep (and whether to keep cellId column).        
    if ~exist('cellId', 'var') || isempty(cellId)      % retrieve spikes for *all* cells (& keep cellId column)        
        spikes = [spikeTimes, spkCellIds, fet];        

    elseif (cellId == 100)      % retrieve spikes for *all* cells (& don't keep cellId column)        
        spikes = [spikeTimes, fet];                
        
    else  % retrieve spikes for 1 specific cell (so can remove cell_id column)
        spk_idx = spkCellIds == cellId;  
        spikes = [spikeTimes(spk_idx), fet(spk_idx,:)];                
    end
    
    if ~isempty(outputType) && ~strcmp(outputType, 'tick')
        Did = dbLookup('Did',  'Gid', Gid);
        spikes(:,1) = dbConvertTimeMeasures(Did, spikes(:,1), 'tick', outputType);
    end    
end
