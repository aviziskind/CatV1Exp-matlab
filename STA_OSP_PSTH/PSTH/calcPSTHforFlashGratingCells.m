function [PSTH_bins, PSTH_vals, PSTH_stats] = calcPSTHforFlashGratingCells(Gid, cellId, psthWindow, opts)

    
    if (nargin < 3) || isempty(psthWindow)
        psthWindow = [-300, 200];
    end
    if (nargin == 4) && isfield(opts, 'nStimMax')
        nStimMax = opts.nStimMax;
    else
        nStimMax = 30;         
    end        
    if (nargin == 4) && isfield(opts, 'stim_ordering')
        stim_ordering = opts.stim_ordering;
    else
        stim_ordering = 'mean';
    end    
    if (nargin == 4) && isfield(opts, 'sortingWindow_ms')
        sortingWindow_ms = opts.sortingWindow_ms;
    else
        sortingWindow_ms = [20, 150];
    end      

    histOpts = struct('psthWindow_ms', psthWindow, 'trialGrouping', 'odd/even');
    [PSTH_bins, stimPSTH_vals_trials, meanFiringRate, bckgSamples] = dbGetCellSpkStimHists(Gid, cellId, histOpts);
    stimPSTH_vals = mean(stimPSTH_vals_trials, 3);
    
    idx_relevant = (PSTH_bins > sortingWindow_ms(1)) & (PSTH_bins < sortingWindow_ms(2)); 
    
    switch stim_ordering
        case 'mean', v = mean(stimPSTH_vals(idx_relevant,:), 1);
        case 'max',  v = max(stimPSTH_vals(idx_relevant,:), [], 1);
    end
    new_order = ord(v, 'descend');
    stimPSTH_vals_ordered = stimPSTH_vals( :, new_order(1:nStimMax) );
    PSTH_vals = mean(stimPSTH_vals_ordered, 2);
    bckgRate = mean(bckgSamples); %[mean(bckgSamples), std(bckgSamples)];
    PSTH_stats = struct('meanRate', meanFiringRate, 'bckgRate', bckgRate );
    
end

%%%  (already do this in dbGetCellSpkStimHists, so removing from here)
%    if cph_separateBeforeAfterTrials && (size(stimPSTH_vals_oe, 3) == 4)
%         bins_before =  ( PSTH_bins <  70 );
%         bins_after  =  ( PSTH_bins >= 70 );                
%         stimPSTH_vals = [mean(stimPSTH_vals_oe(bins_before, :, [1 2]), 3);
%                          mean(stimPSTH_vals_oe(bins_after,  :, [3 4]), 3) ];
%     else            
%         stimPSTH_vals = mean(stimPSTH_vals_oe, 3);
%     end







%
%     switch selectionMethod
%         case 'top',                    
%             % Take top N most preferred stimuli, 
%             N = 20; 
% 
%             new_order = ord(max(stimPSTH_vals, [], 1), 'descend');        
% %             new_order = ord(mean(stimPSTH_vals, 1).*max(stimPSTH_vals, [], 1), 'descend');        
%             topN_PSTHs   =   stimPSTH_vals(:,new_order(1:N));
% %         topN_PSTHs_decor = dcPSTHs(:,new_order(1:N));
%             PSTH_vals = mean( topN_PSTHs, 2 );
% 
% %             osp_mean = reshape(mean(stimPSTH_vals,1), [36, 10, 8]);
% %             osp_max = reshape(max(stimPSTH_vals,[], 1), [36, 10, 8]);
%         case 'repro',
%             
% %             W = .6;
% %             useNStim = false;
% %             psth_data = getPSTHwindowData(Gid, cellId);
% %             [L_bin, R_bin] = getBestPsthWindow(psth_data, W);
% %             bestWindowBins = [L_bin, R_bin];         
%                             
%     end
    
%            meanFiringRate, bckgRate, 
%                 bestWindowBins

%     3;
%     if showWorking
%     figure(50);
%         plotPSTHseries({cumPSTHs(:,1:30)}, '3D');
% %             plotPSTHseries({topN_PSTHs, topN_PSTHs_decor}, '3D');
% %             figure(524); plotThisPSTH(PSTH_bins, PSTH_vals);
%     end
    
% end






%{
%%%% call fig_plotAllStimHists %%%%% 
    getfrm = @getMovieStimulusFrame;
    getfrm('load', 'Gid', Gid);

    [frameStimIds] = getStimulusFrameSequence(Gid, 'OSP');

figure(2);
idx = [new_order(1:20 )];
% idx = [new_order(1801:1805)];
m = 4; n = 5;
hs = zeros(m,n);
ind = 1;
bEdge = binCent2edge(PSTH_bins);
for i = 1:m
    for j = 1:n
%         subplot(m,n,(i-1)*n+j)
        subplot(m,n,ind)
        hs(i,j) = bar(PSTH_bins, stimPSTH_vals(:,idx(ind)),1);
        set(hs(i,j), 'edgecolor', 'b', 'facecolor', 'b');
        ind = ind + 1;
        xlim([bEdge(1), bEdge(end)])
        set(gca, 'xtick', [], 'units', 'pixels')
        p = get(gca, 'position');
        ur = [p(1)+p(3), p(2)+p(4)];
        ll = ur - 16;
        hStim = axes('position', [ll, ur]);
        frmId = find(frameStimIds == idx(ind),1);
        imagesc(getfrm(frmId)); set(hStim , 'xtick', [], 'ytick', []); colormap('gray');
        
%         ylim([0 600])
    end
end
matchAxes('Y', hs(:))







%}


%{
        
            exploreNStimDim = false;
            W = .5;
            minWindowWidth = 2;
            binStart = find(PSTH_bins > 15,1, 'first');
            binEnd   = find(PSTH_bins < 120,1, 'last');
            nStimTh = .75;
            nBins = length(PSTH_bins);
            
            [uori, usp, uph] = dbGetUniqueOriSpPh('Gid', Gid);
            [nOri, nSp, nPh] = deal(length(uori), length(usp), length(uph));
            
            stimPSTH_vals = mean(stimPSTH_vals_oe, 3);
                new_order = ord(mean(stimPSTH_vals, 1), 'descend');        
            stimPSTH_vals_ordered = stimPSTH_vals(:,new_order);
            stimPSTH_vals_cum = cummean(stimPSTH_vals_ordered,2);                
            
            stimPSTH_sums = sum( stimPSTH_vals_ordered, 1);
            nStimMax = find(stimPSTH_sums / stimPSTH_sums(1) <= nStimTh, 1);                    
                        
            if exploreNStimDim
                pvals  = zeros(nBins, nBins, nStimMax);
                slopes = zeros(nBins, nBins, nStimMax);
                psth_areaRatios = zeros(nBins, nBins, nStimMax);

                for nstim = 1:nStimMax
                    cur_psth = stimPSTH_vals_cum(:,nstim);
                    psth_total = sum(cur_psth);

                    for l_bin = binStart : binEnd-1
                        for r_bin = l_bin + (minWindowWidth-1) : binEnd;

                            [pvals(r_bin, l_bin, nstim), slopes(r_bin, l_bin, nstim)] = ...
                                getOspRepVsPsthBinning(stimPSTH_vals_oe, l_bin, r_bin, [nOri, nSp, nPh], cur_psth);

                            frac_area = sum(cur_psth(l_bin:r_bin))/(psth_total);
                            frac_bins = (r_bin-l_bin+1)/nBins;
                            psth_areaRatios(r_bin, l_bin, nstim) = frac_area/(frac_bins);                        
                        end
                    end            
                end
            
                % normalize
                pvals_n = slopes / max(slopes(:));
                psth_area_n = psth_areaRatios / max(psth_areaRatios(:));
                T = (W*pvals_n + (1-W)*psth_area_n);

                [max_T, idx_maxT] = maxElement(T);
                [R_bin_best, L_bin_best, nStim_best] = dealV( idx_maxT );

                bestWindow = [L_bin_best, R_bin_best];
                PSTH_vals = stimPSTH_vals_cum(:,nStim_best);            
            else

                pvals  = zeros(nBins, nBins);
                slopes = zeros(nBins, nBins);
                psth_areaRatios = zeros(nBins, nBins);

                for l_bin = binStart : binEnd-1
                    for r_bin = l_bin + (minWindowWidth-1) : binEnd;

                        [pvals(r_bin, l_bin), slopes(r_bin, l_bin, nstim)] = ...
                            getOspRepVsPsthBinning(stimPSTH_vals_oe, l_bin, r_bin, [nOri, nSp, nPh]);

                        frac_area = sum(cur_psth(l_bin:r_bin))/(psth_total);
                        frac_bins = (r_bin-l_bin+1)/nBins;
                        psth_areaRatios(r_bin, l_bin, nstim) = frac_area/(frac_bins);                        
                    end
                end                            
            
                % normalize
                pvals_n = slopes / max(slopes(:));
                psth_area_n = psth_areaRatios / max(psth_areaRatios(:));
                T = (W*pvals_n + (1-W)*psth_area_n);

                [max_T, idx_maxT] = maxElement(T);
                [R_bin_best, L_bin_best, nStim_best] = dealV( idx_maxT );

                bestWindow = [L_bin_best, R_bin_best];
                PSTH_vals = stimPSTH_vals_cum(:,nStim_best);                     
                
                
                
                
                
            end
%}        