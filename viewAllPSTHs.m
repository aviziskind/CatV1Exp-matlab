%{
load indivCells_movie_fg


%}


% inds = findInStructArray(movieGroups_fg, 'frameLength_ms', [], @(x) x < 20);
% for i = 1:length(inds)
%     Gid = movieGroups_fg(inds(i)).Gid;
%     cellIds = movieGroups_fg(inds(i)).cellIds;
%     for c = cellIds
%         varname = getName('celldata', Gid, c);
%         v = eval(varname);
%         plotThisPSTH(v.PSTH); 
%         fprintf('Group %d, cell %d\n', Gid, c);
%         3;
%     end    
% end



allFunnyGids = [ 2625 2665 2677 2723 2737 2747 2775 2791 2805 2887 2921 2957 2981 3047 3061 3077 3091];
inds = findInStructArray(movieGroups_fg, 'Gid', [], @(x) any(x == allFunnyGids));
for i = 1:length(inds)
    Gid = movieGroups_fg(inds(i)).Gid;
    cellIds = movieGroups_fg(inds(i)).cellIds;
    for c = cellIds
        varname = getName('celldata', Gid, c);
        v = eval(varname);
        plotCellData(v); 
        fprintf('Group %d, cell %d\n', Gid, c);
        3;
    end    
end
