function fig_measures

    function [r1, r2] = getTCs(GC)

        [Gid1, cellId1, Gid2, cellId2, ori_i, sp_i] = elements(GC);
        S = load('flashedGratingCells_DB_all');
        allCells = S.allCells;
        gids = [allCells.Gid];
        cellids = [allCells.cellId];
        id1 = find(gids == Gid1 & cellids == cellId1, 1);
        id2 = find(gids == Gid2 & cellids == cellId2, 1);
        OSP1 = allCells(id1).R;
        OSP2 = allCells(id2).R;
        r1 = squeeze(OSP1(ori_i, sp_i,:));
        r2 = squeeze(OSP2(ori_i, sp_i,:))*1.5;    
        r1 = double(r1);
        r2 = double(r2);
        
    end

    wrp = @(x) x([1:end, 1]);
    shftR = @(x) x([end, 1:end-1]);
    shftRn = @(x,n) x([end-n+1:end, 1:end-n]);

m_cc = 1;
m_rho = 2;
m_dphi = 3;
m_dF1 = 4;

doFigs = [m_dF1];

phs = 0:45:315;
nPh = length(phs);
phs_360 = 0:45:360;

if any(doFigs == m_cc)

    GC = [1875 1 1877 0 3 4];
    
%     GC = [2272           5        2272           6           7           1];
    [r1, r2] = getTCs(GC);

    r2 = r2*2;
    figure(1);
    h1 = plot(phs_360, wrp(r1), 'bo-', phs_360, wrp(r2), 'gs-');
    
    r1ms = r1 - mean(r1);
    r2ms = r2 - mean(r2);    
    xlim([0, 360]);
    set(gca, 'xtick', phs_360);
    xlabel('\phi');
    ylabel('Firing rate');
    title('Phase tuning curves');    
    drawHorizontalLine(0, 'linestyle', ':')
    
    figure(2);
    h2 = plot(phs_360, wrp(r1ms), 'bo-', phs_360, wrp(r2ms), 'gs-');
    xlim([0, 360]);
    set(gca, 'xtick', phs_360);
    xlabel('\phi');
    ylabel('Firing rate above mean');
    title('Phase tuning curves (mean subtracted)');    
    drawHorizontalLine(0, 'linestyle', ':')
    
    matchAxes('Y', [h1 h2])
    cc = pearsonR(r1, r2)                   %#ok<*NASGU>
    dphi = deltaPhi(phs, r1, r2, 'angle')
    
    
end


if any(doFigs == m_rho)

    GC2 = [2903, 0, 4774, 3, 14, 7];  % pathological case: high cc, low rho.    
    GC3 = [2975, 6, 4746 1, 21, 2];% path case 2: high rho, low cc;
    
    GC = [1875 0 1875  1 6 4];
    
    [r1, r2] = getTCs(GC);

    r2 = r2*2;

    R1 = tiedrank(r1);
    R2 = tiedrank(r2);
    
    R1ms = R1 - mean(R1);
    R2ms = R2 - mean(R2);    
    
    figure(1);
    h1 = plot(phs_360, wrp(r1), 'bo-', phs_360, wrp(r2), 'gs-');    
    xlim([0, 360]);
    set(gca, 'xtick', phs_360);
    xlabel('\phi');
    ylabel('Firing rate');
    title('Phase tuning curves');    
    drawHorizontalLine(0, 'linestyle', ':')
    
    figure(2);
    h2 = plot(phs_360, wrp(R1), 'bo-', phs_360, wrp(R2), 'gs-');
    xlim([0, 360]);
    set(gca, 'xtick', phs_360);
    xlabel('\phi');
    ylabel('Rank');
    title('Phase tuning curves (ranked)');    
    drawHorizontalLine(0, 'linestyle', ':')

    figure(3);
    h3 = plot(phs_360, wrp(R1ms), 'bo-', phs_360, wrp(R2ms), 'gs-');
    xlim([0, 360]);
    set(gca, 'xtick', phs_360);
    xlabel('\phi');
    ylabel('Rank - mean rank');
    title('Phase tuning curves (ranked)');    
    drawHorizontalLine(0, 'linestyle', ':')
    
    
    matchAxes('Y', [h2 h3])
%     cc = pearsonR(r1, r2)
%     dphi = deltaPhi(phs, r1, r2, 'angle')
    rho = spearmanRho(r1, r2) %#ok<*NOPRT>
    
    
end



if any(doFigs == m_dphi)
    MODE = 'cross-correlation';
    % MODE = 'angle';

    GC = [1875 0 1875 1  6 4];
    [r1, r2] = getTCs(GC);

    ccs = zeros(1,nPh);
    r2_0 = r2;
    for i = 1:nPh    
        ccs(i) = pearsonR(double(r1), double(r2));
        r2 = shftR(r2);
    end
    ccs = wrp(ccs);
    r2 = r2_0;

    figure(5); clf; 
    subplot(1,2,1);
    h1 = plot(phs_360, wrp(r1), 'bo-'); hold on;
    h2s = plot(phs_360, wrp(r2), 'g.:');
    h2 = plot(phs_360, wrp(r2), 'gs-');
    xlim([0, 360]);

    set(gca, 'xtick', phs_360);
    xlabel('\phi');
    ylabel('Firing rate');
    title('Phase tuning curves');

    subplot(1,2,2);
    hccs = plot(phs_360, ccs, 'k.:'); hold on
    hcc1 =plot(phs(1), ccs(1), 'o', 'color', 'k');
    cc_tick = [0:45:180, 135:-45:0];
    set(gca, 'xtick', phs_360);
    xlim([0, 360]);
    title('Correlation coefficient');

    axis(axis);
    xlabel('\Delta\phi');
    ylabel('r');
    set(gca, 'xticklabel', cc_tick);

    r2 = r2_0;
    for i = 1:nPh
        r2 = shftR(r2);
        set(h2, 'ydata', wrp(r2));
        cc = pearsonR(double(r1), double(r2));
        set(hcc1, 'xdata', phs_360(i+1), 'ydata', ccs(i+1));
            3;
    end

    deltaPhi(phs, r1, r2, MODE);

end
    
if any(doFigs == m_dF1)
    
    
    GC = [2302 1 2302  4  22 4];
    [r1, r2] = getTCs(GC);
    r2 = r2*2;
           
    r2 = shftRn(r2,3);
%     [phi1] = getF1phase(deg2rad(phs), r1, 2*pi);
%     [phi2] = getF1phase(deg2rad(phs), r2, 2*pi);
    [r1_F1, t1, f_cos1] = firstNHarmonics(r1, 1, 100, [], phs);
    [r2_F1, t2, f_cos2] = firstNHarmonics(r2, 1, 100, [], phs);
    
%     [phi1, f_cos1, t1] = getF1phase(deg2rad(phs), r1, 2*pi);
%     [phi2, f_cos2, t2] = getF1phase(deg2rad(phs), r2, 2*pi);
%     ph1 = rad2deg(2*pi-phi1);
%     ph2 = rad2deg(2*pi-phi2);
    
    
    figure(1); clf;
    h1 = plot(phs_360, wrp(r1), 'bs-'); hold on;
    h2 = plot(phs_360, wrp(r2), 'gs-');
    xlim([0, 360]);
    ylims = ylim;
    ylim([-5, ylims(2)]);
    set(gca, 'xtick', 0:90:360)
    xlabel('\phi');
    ylabel('firing rate');

    3;
    plot(t1, f_cos1, 'b:');
    plot(t2, f_cos2, 'g:');
    plot(phs_360, wrp(r1_F1), 'bo', 'markerfacecolor', 'b');
    plot(phs_360, wrp(r2_F1), 'go', 'markerfacecolor', 'g');
    
    3;
    [mx1, ind1] = max(r1_F1);
    [mx2, ind2] = max(r2_F1);
    phi1 = phs(ind1);
    phi2 = phs(ind2);
%     drawVerticalLine([phi1, phi2], 'LineStyle', ':');
%     y = mean([mx1, mx2]);
    y = mx1*1.1;%mean([mx1, mx2]);
    [fx, fy] = dsxy2figxy(gca, [phi1, phi2], [y y]);
    annotation('doublearrow', fx, fy);    
    
    phi1_rad = deg2rad(phi1);
    phi2_rad = deg2rad(phi2);
    figure(2); clf;
    h1p = polar(deg2rad(phs_360), wrp(r1)', 'bo-'); hold on;
    h2p = polar(deg2rad(phs_360), wrp(r2)', 'gs-');
    set([h1p, h2p], 'linewidth', 1);
    R1 = max(r1)*.8;
    R2 = max(r2)*.8;
    warning('off', 'MATLAB:plot:DeprecatedV6Argument')
    h1q = quiver('v6', 0,0, R1*cos(phi1_rad), R1*sin(phi1_rad), 'b'); %#ok<*NOV6>
    h2q = quiver('v6', 0,0, R2*cos(phi2_rad), R2*sin(phi2_rad), 'g');
    set([h1q h2q], 'linewidth', 2.5)
    3;
    
%     deltaPhi(phs, r1, r2, MODE)

end


end


