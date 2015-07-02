% function plotMeasuresVsMeasures

    cc = 'cc';
    rho = 'rho';
    tau = 'tau';
    dphi = 'dphi';
    dF1 = 'dF1';
%{
    fg_cmpDatafile = [CatV1Path 'flashedGratingComparisonData.mat'];     
    S_fg = load(fg_cmpDatafile); 
    [fg_pairData, S_fg, pairTypes, measures, locations] = deal(S_fg.pairData, S_fg.allStatsC, S_fg.pairTypes, S_fg.measures, S_fg.locations);
    dg_cmpDatafile = [CatV1Path 'driftingGratingComparisonData.mat'];     
    S_dg = load(dg_cmpDatafile); 
    [dg_pairData, S_dg, pairTypes, measures, locations] = deal(S_dg.pairData, S_dg.allStatsC, S_dg.pairTypes, S_dg.measures, S_dg.locations);

    %}
%     
S = S_dg;
pairData = dg_pairData;

allNphs = [pairData.n_phases];
N = 1000;    
Nphs = unique(allNphs);
% Nphs = Nphs(1);
idx = arrayfun(@(n) find(allNphs == n, N, 'first'), Nphs, 'un', 0);    
% idx4 = find(allNphs == 8, N, 'first')
% idx8 = find(allNphs == 8, N, 'first');    

% measures = 'cc'    'rho'      'dphi'    'dF1';
measureLabels = {'r', '\rho', '\Delta\phi', '\DeltaF1'};

% doFigs = {rho, cc; dphi, dF1; cc, dF1; rho, tau};
doFigs = {rho, dF1; cc, dF1};
cols = 'bg';

for fig_i = 1:size(doFigs,1);
    m1 = doFigs{fig_i,1};
    m2 = doFigs{fig_i,2};
    figure(fig_i); clf;
    
    m_ind1 = find(strcmp(m1, measures));
    m_ind2 = find(strcmp(m2, measures));
    
    for ni = 1:length(Nphs)
        m1_ph{ni} = S{1,m_ind1}.val(idx{ni}); %#ok<*AGROW>
        m2_ph{ni} = S{1,m_ind2}.val(idx{ni});
    end
%     m1_ph8 = S{1,m_ind1}.val(idx8);
%     m2_ph8 = S{1,m_ind2}.val(idx8);
    m1_lims = S{1,m_ind1}.binEdges([1,end]);
    m2_lims = S{1,m_ind2}.binEdges([1,end]);
    
    m1_lims = iff(m1_lims(2) >= 180, [-1, 181], [-1, 1]);
    m2_lims = iff(m2_lims(2) >= 180, [-1, 181], [-1, 1]);
    
    tick1 = iff(m1_lims(2) >= 180, [0:45:180], [-1:.5:1]);
    tick2 = iff(m2_lims(2) >= 180, [0:45:180], [-1:.5:1]);
    
    slope = iff( xor(m1_lims(2) >= 180,  m2_lims(2) >= 180), -1, 1) ;
    
    for ni = 1:length(Nphs)
        h(ni) = plot(0,0, [cols(ni) 'o'], 'markerfacecolor', cols(ni), 'markersize', 3); hold on;
    end    
    for ni = 1:length(Nphs)
        plot(m1_ph{ni}, m2_ph{ni}, [cols(ni) 'o'], 'markersize', 2);
    end
    set(h(1:length(Nphs)), 'visible', 'off')
    X = [tick1(1) tick1(end)];
    Y = [tick2(1) tick2(end)]; if slope == -1, Y = fliplr(Y); end
    plot(X, Y, 'k')
%     slope
    
    xlabel(measureLabels(m_ind1));
    ylabel(measureLabels(m_ind2));
    set(gca, 'xtick', tick1, 'ytick', tick2);
    daspect([1, diff(m2_lims)/diff(m1_lims), 1]);
    axis([m1_lims m2_lims])
    if length(Nphs) == 1
        title(['n = ', num2str(Nphs)]);
    else
        hL = legend(legendarray('n = ', Nphs), 'location', 'bestoutside');
        set(hL, 'fontsize', 9);
    end
    
%     axis equal tight;
end    
    
    

% end




% cmp_gids = fg_pairData(ind).Gids;
% cmp_cellids = fg_pairData(ind).cellIds;
% 
% 
% id1 = find(gids == cmp_gids(1) & cmp_cellids(1) == 0);
% id2 = find(gids == cmp_gids(2) & cmp_cellids(2) == 3);
% R1 = allOSPs(id1).R;
% R2 = allOSPs(id1).R;
% mR1 = mean(R1, 3);  fR1 = mR1/max(mR1(:));
% mR2 = mean(R2, 3);  fR2 = mR2/max(mR2(:));
% [tmp, oriSpInds] = maxElement(fR1 .* fR2);
% ph1 = double(squeeze(R1(oriSpInds(1),oriSpInds(2),:)));
% ph2 = double(squeeze(R2(oriSpInds(1),oriSpInds(2),:)));
% phs = 0:45:315;
% pearsonR(ph1, ph2);
% r1 = tiedrank(ph1); r2 = tiedrank(ph2);
% 
% pearsonR(double(ph1), double(ph2));