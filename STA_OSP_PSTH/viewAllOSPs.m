function viewAllOSPs(allOSPs)

    stimType = 'flashed'; % flash or drifting.

    if nargin < 1
        S = load([ stimType 'GratingOSPs.mat']);
        allOSPs = S.allOSPs;
    end

%     sortMode = 'manualRank';
    sortMode = 'stats';
    
    switch sortMode
        case 'manualRank'
            rnks = [allOSPs.manualRank];            
%             inds = ord(rnks);
            inds = find(rnks == 3);
            allOSPs = allOSPs(inds);
            ws = rnks(inds);
            
        case 'autoRank'
            
            
        case 'stats'
            % select by >10 spfs.            
            nsp = arrayfun(@(x) length(x.sp), allOSPs);
            inds = (nsp >= 10);
            allOSPs = allOSPs(inds);

            % select by all 3 pvalues < .01
            allstats = [allOSPs.stats];

%             fld = 'orientationSelectivePval';
%             fld = 'orientationReproduciblePval';
%             fld = 'spatFreqReproduciblePval';
%             fld = 'spatFreqReproduciblePval';
%             pval = [allstats.fld];
%             ori_sel_pval = [allstats.orientationSelectivePval];
%             ori_rep_pval = [allstats.orientationReproduciblePval];
%             spf_rep_pval = [allstats.spatFreqReproduciblePval];
            resp_pval = -log10([allstats.responseReproduciblePval]);
            mx = max(resp_pval(resp_pval<inf));
            resp_pval(resp_pval==inf)=mx;
            
%             scores = ori_sel_pval+ori_rep_pval+spf_rep_pval;
%             inds = (ori_sel_pval < .02) & (ori_rep_pval < .02) & (spf_rep_pval < .02);

            inds = ord(resp_pval, 'ascend');
            allOSPs = allOSPs(inds);
            
            ws = resp_pval(inds);
            % sort by max firing rate.
%             max_r = arrayfun(@(x) max(x.OSP(:)), allOSPs);            
%             [scores_sorted, inds] = sort(max_r, 'descend');
%             allOSPs = allOSPs(inds);
            
        case 'spf_pval'
            
    end
    
    m = 4; n = 5;
    gridSubPlot(m,n,  [1 20]);
    for osp_i = 1:length(allOSPs)        
%         subplot(m,n, mod(osp_i-1, m*n)+1);
        lastOne = gridSubPlot; 
       
        S = allOSPs(osp_i);
        [groupId, cellId] = deal( S.GroupId, S.cellId );
%         spfTuningCurve = S.spfTuningCurve;
        imageOSP(S, 'mean:ph', 'SOP', 'noTicks', 'noLabels'); %  colorbar;

%         nullSpfCurve = ones(size(spfTuningCurve)) * sum(spfTuningCurve)/length(spfTuningCurve);
%         cla;
%         bar(spfTuningCurve); hold on;
%         xlim(.5+[0, length(spfTuningCurve)])
%         drawHorizontalLine( nullSpfCurve(1), 'color', 'r');
%         [spf_sel,spf_pval, cstat] = histChiSqrTest(spfTuningCurve, nullSpfCurve, alpha);
        
        axis normal
        set(gca, 'xtick', [], 'ytick', []); 
                
%         w = findHowWellTunedOSP(OSP);
        w = ws(osp_i);% ws(osp_i);%stats.(fld);
%         w = scores_sorted(osp_i);
        title(sprintf( '(%d/%d) %d, %d. w = %3.2g', osp_i, length(allOSPs), groupId, cellId, w));
        set(gca, 'xtick', []);
        3;
        if mod(osp_i, m*n) == 0
            input('Press <enter> to continue');
        end    
%         if lastOne, break, end;        
    end

end

%  some nice (fg) cells
% 4538, 0
% 2771, 9
% 4482, 0
% 4494, 0

% some nice (dg) cells
% 1145, 0 *
% 1140, 0
% 1171, 0
% 4047, 0
% 1131, 4 (very sparse)

%{
    if redo 
    
        if nargin < 1
            S = load([ stimType 'GratingPairs.mat']);
            allOSPs = S.allOSPs;
        end
        excludeNonMovieStims = false;
        if excludeNonMovieStims
            inds = findInStructArray(allOSPs, 'ori', [], @(x) length(x) > 30);
            allOSPs = allOSPs(inds);
        end
        alpha = .05;
        % add sp_pval field
        for i = 1:length(allOSPs);
            R = allOSPs(i).R;
            [tmp, indmax] = maxElement(R);
            [ori_peak_ind, sp_peak_ind] = elements(indmax);
            spfTuningCurve = mean( R(ori_peak_ind, :,:) ,3); % at peak ori, sum across phases
            spfTuningCurve = sum(spfTuningCurve(ind),1);
%             spfTuningCurve = spfTuningCurve/sum(spfTuningCurve);
%             
            nullSpfCurve = ones(size(spfTuningCurve)) * sum(spfTuningCurve)/length(spfTuningCurve);
            allOSPs(i).spfTuningCurve = spfTuningCurve;
            if (length(nonnans(spfTuningCurve)) > 1)        
                [spf_sel,spf_pval] = histChiSqrTest(spfTuningCurve, nullSpfCurve, alpha);
                allOSPs(i).spf_pval = spf_pval; 
            end        
        end
    
        save('allOSPs_tmp.mat', 'allOSPs');
    
    else        
        S = load('allOSPs_tmp.mat');
        allOSPs = S.allOSPs;
    end
%}