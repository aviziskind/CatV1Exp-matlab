
for i = 16:50
    grp = movieGroups( bandNoiseInds(i) );
    Gid = grp.Gid;
    cells = grp.cellIds;
    
    figure(i); clf;
    for cell_i = 1:length(cells)
        nm = ['celldata__Group_' num2str(Gid) '__Cell_' num2str(cells(cell_i)) ];
        data = eval(nm);
        subplot(2, length(cells), cell_i);
        if ~isstruct(data), continue; end
        imagesc(data.STAs.STA(:,:,1)); set(gca, 'xtick', [], 'ytick', []);
        axis equal tight;

        subplot(2, length(cells), length(cells)+cell_i);
        imagesc(data.STAs.STA(:,:,2)); set(gca, 'xtick', [], 'ytick', []);
        axis equal tight;
        
    end
    
    input('');
    
end