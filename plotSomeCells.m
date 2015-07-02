Gids = [ 2625 2665 2677 2723 2737 2747 2775 2791 2805 2887 2921 2957 2981 3047 3061 3077 3091];

for Gi = 1:length(Gids)
    grp = movieGroups_fg(findInStructArray(movieGroups_fg, 'Gid', Gids(Gi)));
    cellIds = grp.cellIds;
    
    for ci = 1:length(cellIds)
        varname = getName('celldata', Gids(Gi), cellIds(ci));
        v = eval(varname);

        figure(1);
        plotCellData(v);  axis square;
        3;
    end
    
end