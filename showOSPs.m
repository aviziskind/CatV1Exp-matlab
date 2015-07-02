
% show multiple OSP profiles
inds1 = findInStructArray(movieGroups_fg, 'cellIds', [], @(x) length(x) >1);
inds2 = findInStructArray(movieGroups_fg, 'spPh_deg', [], @(x) length(x) >5);
inds = intersect(inds1, inds2);
for i = 1:length(inds)
    grp = movieGroups_fg( inds(i) );
    Gid = grp.Gid;
    cells = grp.cellIds;
    
    figure(1); clf;
    for cell_i = 1:length(cells)
        nm = getName('celldata', Gid, cells(cell_i));
        data = eval(nm);
        h(cell_i) = subplot(1, length(cells), cell_i);
        if ~isstruct(data), continue; end
        imageOSP(data.OSP.R, 'mean:Phase'); colorbar;
        
        xlabel(num2str(cells(cell_i)));
    end 
%     matchAxes('C', h(1:length(cells)));
    
    suptitle_2([' group ' num2str(Gid)] );
    input('');
    
end