allCells = [1857 2
1857 3
1931 3
1979 1
2022 1
2022 3
2022 4
2066 1
2066 2
2216 1];

downSampFactor = 2;
pth = 'C:\ExperimentDB\MID\Cells\';
for ci = 1:size(allCells, 1)
    Gid = allCells(ci,1);
    cellId = allCells(ci,2);
    
    mid_fileName_ustim = getName('MID_file', Gid, cellId, downSampFactor, 'uStim');    
    [pth_orig fileName_u, ext] = fileparts(mid_fileName_ustim);
    mid_fileName_ustim = [pth, fileName_u, ext];
    
    mid_fileName_all = getName('MID_file', Gid, cellId, downSampFactor, 'all');
    [pth_orig fileName_a, ext] = fileparts(mid_fileName_all);
    mid_fileName_all = [pth, fileName_a, ext];
    
    S_u = load(mid_fileName_ustim);
    S_a = load(mid_fileName_all);
    
    s = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
    
    figure(ci); clf;
    h_ax_sta = subplot(1,4,1); h(1) = imagesc(s.STAs.STA); title(sprintf('Gid = %d, cellId = %d, STA', Gid, cellId)); axis square equal tight;
    h_ax(1) = subplot(1,4,2); h(1) = imagesc(S_u.MID); title(sprintf('Gid = %d, cellId = %d, uStim', Gid, cellId)); axis square equal tight;
    h_ax(2) = subplot(1,4,3); h(2) = imagesc(S_a.MID); title(sprintf('Gid = %d, cellId = %d, all', Gid, cellId)); axis square equal tight;
    h_ax(3) = subplot(1,4,4); h(3) = imagesc(S_a.MID-S_u.MID); colorbar; title(sprintf('Gid = %d, cellId = %d, diff', Gid, cellId)); axis square equal tight;
    L_sta = get(h_ax_sta, 'clim'); set(h_ax_sta, 'clim', max(abs(L_sta))*[-1, 1]);
    clims = get(h_ax, 'cLim');
    L = max(abs([clims{:}]));
    set(h_ax, 'clim', [-L, L]);

    3;
end