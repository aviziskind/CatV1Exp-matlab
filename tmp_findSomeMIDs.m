function tmp_findSomeMIDs

Gid = 4502;
cellIds = [1 4 5];
Gids = Gid*ones(1, length(cellIds));
frameModes = {'uStim', 'all'};
opt.downSmpFactor = 4; % degree of downsampling.        

for ci = 1:length(cellIds)
    cellId = cellIds(ci);
    for fi = 1:length(frameModes);
        opt.frameMode = frameModes{fi};        
        
        MID = mid_findMostInfoDim(Gid, cellId, opt);
    end
end

