
allCells = ...
[4476,0; 4470,2; 2544,3; 4878,4;
 4508,2; 4522,3; 2236,0; 2066,1;
 4706,1; 4898,3; 4892,3; 5164,4;
 2244,5; 4476,2; 2434,3; 4936,1; 
 4518,4; 
 4482,0;  
 4474,0;
 4506,0;
 4990,4;
 4898,3;
 3057,0;
 2288,1;];
stimType = 'flashed'; % flash or drifting.

% if nargin < 1
    S = load([ stimType 'GratingCells.mat']);
    allOSPs = S.allOSPs;
% end

    m = 4; n = 4;
    gridSubPlot(m,n,  [1 20]);

gids = [allOSPs.GroupId];
cellids = [allOSPs.cellId]    ;
    
for osp_i = 1:length(allCells)        
%         subplot(m,n, mod(osp_i-1, m*n)+1);
        lastOne = gridSubPlot; 
        set(gcf, 'color', 'w');
        id = find( (gids == allCells(osp_i,1)) & (cellids == allCells(osp_i,2) ));
        
        S = allOSPs(id);
        [groupId, cellId] = deal( S.GroupId, S.cellId );
%         spfTuningCurve = S.spfTuningCurve;
        imageOSP(S, 'mean:ph', 'SOP', 'noTicks', 'noLabels'); %  colorbar;

%         nullSpfCurve = ones(size(spfTuningCurve)) * sum(spfTuningCurve)/length(spfTuningCurve);
%         cla;
%         bar(spfTuningCurve); hold on;
%         xlim(.5+[0, length(spfTuningCurve)])
%         drawHorizontalLine( nullSpfCurve(1), 'color', 'r');
%         [spf_pval, spf_sel, cstat] = histChiSqrTest(spfTuningCurve, nullSpfCurve, alpha);
        
        axis normal
        set(gca, 'xtick', [], 'ytick', []); 
                
%         w = findHowWellTunedOSP(OSP);
%         w = ws(osp_i);% ws(osp_i);%stats.(fld);
%         w = scores_sorted(osp_i);
        title(sprintf('%d : %d', groupId, cellId));
        set(gca, 'xtick', []);
        3;
        if mod(osp_i, m*n) == 0
            input('Press <enter> to continue');
        end    
%         if lastOne, break, end;        
end