function stats = calcStatsFromOSPfull(R, R_full, f, Gid, bckgRate)

    repTest = 'spearman'; % options: 'spearman', 'regression'.

    if iscell(R)
        stats = cellfun(@(R1, R_full1, f1, uori1) ...
                calcStatsFromOSPfull(R1, R_full1, f1, uori1, Gid, bckgRate), R, R_full, f, uori, 'un', 0);
        return;
    end

    global showWorking;
    if isempty(showWorking)
        showWorking = false;
    end        
        
    if ~exist('f', 'var') || isempty(f)
        f = 1;
    end    
    if strcmp(class(R_full), 'uint8')                
        if ~isstruct(f), f = struct('scale', f, 'idxnan', []); end
        R_full = single(R_full)/f.scale;
        R_full(f.idxnan) = nan;
    end
    stats = struct;
    
    if exist('Gid', 'var') && ~isempty(Gid)
        gratingType = flashedOrDrifting(Gid, 's');        
    else
        gratingType = [];
    end
    [nOri, nSpf, nPh, nTrialsMax] = size(R_full); 
    
    % REPRODUCITBILITY TESTS
    % divide the data in half, and check that orientation tuning curves
    % from each half are significantly correlated.

    %     phase_action = {'average'};  %options: {'average', 'itemize', 'max'}
    trial_action = 'odd/even'; %options: 'even/odd', 'first/second'

    showBckgInReproPlots = true;
%     if exist('bckgRate', 'var') && (length(bckgRate) > 1) && ~isnan(bckgRate(1)) 
%         [bckgMean, bckgStd] = elements(bckgRate);

    
    %  TEST RESPONSE REPRODUCITBILITY 
%     [oriSpf_rep, oriSpf_rep_pval, oriSpf_rep_str] = deal(0, 1, 1, 0);
%     [oriSpf_rep_top, oriSpf_rep_pval_top, oriSpf_rep_str_top] = deal(0, 1, 1, 0);
%     [oriSpf_rep_smth, oriSpf_rep_pval_smth, oriSpf_rep_str_smth] = deal(0, 1, 1, 0);
        
        
    inclCriteria = {'all_over', @ge, 0}; % options: 'all',  'all_over' (ie. all stimuli where response was over x % of max); 'top' (ie. top x % of stimuli)
%     inclCriteria = {'all'}; 
    doRepTests = {'ori_sp_avPhase', 'ori_sp_maxPhase', 'ori_sp_ph', 'ori_sp_smoothed'};

        
    %1. test ori-sp response reproducibility (averaging over phases)    
    testType = 'ori_sp_avPhase';
    if any(strcmp(testType, doRepTests))        
        [showWorkingHere, figId] = testIfShowWorkHere(showWorking, testType, 501);        
        [rep, rep_pval, rep_str] = ...
            doReproducibilityTest(repTest, R_full, 'average', trial_action, inclCriteria, showWorkingHere, figId, testType);
        stats.(['rep_' testType '_pval']) = rep_pval;
        stats.(['rep_' testType '_str']) = rep_str;
    end
                    
    %2. test ori-sp response reproducibility (taking max of phases)
    testType = 'ori_sp_maxPhase';
    if any(strcmp(testType, doRepTests))
        [showWorkingHere, figId] = testIfShowWorkHere(showWorking, testType, 502);

        [rep, rep_pval, rep_str] = ...
            doReproducibilityTest(repTest, R_full, 'max', trial_action, inclCriteria, showWorkingHere, figId, testType);        
        stats.(['rep_' testType '_pval']) = rep_pval;
        stats.(['rep_' testType '_str']) = rep_str;        
    end    

    %3. test ori-sp-ph response reproducibility
    testType = 'ori_sp_ph';
    if any(strcmp(testType, doRepTests))
        [showWorkingHere, figId] = testIfShowWorkHere(showWorking, testType, 503);
        
        [rep, rep_pval, rep_str] = ...
            doReproducibilityTest(repTest, R_full, 'itemize', trial_action, inclCriteria, showWorkingHere, figId, testType);        

        stats.(['rep_' testType '_pval']) = rep_pval;
        stats.(['rep_' testType '_str']) = rep_str;                    
    end
    
    %4. if have not smoothed over phases, and is drifting grating, also
    %compute smoothed reproducibility (compressing over phases).
    testType = 'ori_sp_ph_smoothed';    
    if any(strcmp(testType, doRepTests)) && strcmp(gratingType, 'drifting') && (nPh > 15)            
        R_full_smoothed = compressOSP_Phs(R_full);        
        [showWorkingHere, figId] = testIfShowWorkHere(showWorking, testType, 504);

        [rep, rep_pval, rep_str] = ...
            doReproducibilityTest(repTest, R_full_smoothed, 'average', trial_action, inclCriteria, showWorkingHere, figId, testType);        

        stats.(['rep_' testType '_pval']) = rep_pval;
        stats.(['rep_' testType '_str']) = rep_str;        
    end
    
               
    
    % 4. test RESPONSE SIZE SIGNIFICANCE
    % test that the response at the preferred orientation/sp is >3 * std.dev +
    % spontaneous firing rate.    
    [r_size, r_size_frac] = deal([]);
    showRespSize = (islogical(showWorking) && showWorking) || (ischar(showWorking) && ~isempty(strfind(showWorking, 'resp_size')));
    if 0 && exist('bckgRate', 'var') && (length(bckgRate) > 1) && ~isnan(bckgRate(1)) 
        R_ori_spf = mean(R,3);
        [tmp, indmax] = maxElement(R_ori_spf);
        [ori_peak_ind, spf_peak_ind] = elements(indmax);            
                
        [bckgMean, bckgStd] = elements(bckgRate);
        respToBestOriSp = mean( R(ori_peak_ind, spf_peak_ind, :) );
        r_size_frac = (respToBestOriSp-bckgMean) / bckgStd;
        r_size = r_size_frac > 3;
%             stats.responseSize = (respToBestOri > bckgMean + 3*bckgStd);
%             stats.responseSizeFrac = (respToBestOri-bckgMean) / bckgStd;

        if showRespSize
            figure(112); clf;
            line([0;1], [1;1]*respToBestOriSp, 'color', 'b')
            line([0;1], [1;1]*bckgMean, 'color', 'r')
            line([0;1], [1;1]*bckgMean+3*bckgStd, 'color', 'r', 'linestyle', ':')
            hold on; plot(0,0);
            title( sprintf('Significant response : h = %d. frac = %2.4g', r_size, r_size_frac));
        end            
    end
    stats.responseSize = r_size;
    stats.responseSizeFrac = r_size_frac;

    
    
end
    



%----------------------------------------------------------------------
function [showWorkingHere, figId] = testIfShowWorkHere(showWorking, testType, figId)
    showWorkingHere = (islogical(showWorking) && showWorking) || (ischar(showWorking) && ~isempty(strfind(showWorking, testType)));
    if (ischar(showWorking) && ~isempty(strfind(showWorking, testType))) 
        s = strtok(showWorking, testType);
        id = str2double(s);
        if ~isnan(id)
            figId = id;
        end
    end
end
    
%----------------------------------------------------------------------
function [rep, rep_pval, rep_str] = doReproducibilityTest(repTestType, Rf, ph_action, tr_action, inclCriteria, showTF, fig_id, lbl)
    if showTF
        figure(fig_id); clf;
    end
    alpha = .01;
    [oriSpfTuning_1, oriSpfTuning_2] = rearrangeTrials(Rf, inclCriteria, ph_action, tr_action);
    switch repTestType
        case 'spearman',
            [rho, rep_pval, rep_str] = myCorr(oriSpfTuning_1(:), oriSpfTuning_2(:), 'type', 'spearman', 'tail', 'gt');
            rep = rep_pval < alpha;
        case 'regression',     
            [rep, rep_pval, rep_str] = signRegressionTest(oriSpfTuning_1(:), oriSpfTuning_2(:), alpha, iff(showTF, 1, []));
    end
    if showTF && exist('lbl', 'var') && ~isempty(lbl);
%             xlabel('odd trials'); ylabel('even trials');
        xlabel(lbl, 'interpreter', 'none');
    end

end


%----------------------------------------------------------------------
function [rep, rep_pval, rep_str] = signRegressionTest(x, y, alpha, showWorkingFlag)

    [rep, rep_pval, rep_slope] = regressionSlopeTtest(x, y, alpha, '+', showWorkingFlag);
%     if rep_slope < 0
%         [tmp1, rep_pval2] = regressionSlopeTtest(x, y, alpha, '-');
%         rep_pval_sgn = 1/rep_pval2;
%     elseif rep_slope >= 0
%         rep_pval_sgn = rep_pval;
%     end

    rep_str = rep_slope;
    if rep_slope < 0
        rep_str = max(rep_slope, 1/rep_slope);
    elseif rep_slope > 0
        rep_str = min(rep_slope, 1/rep_slope);
    end                            

%     if ~isempty(showWorkingFlag) && showWorkingFlag
%         xlims = xlim; ylims = ylim; L = max(xlims(2), ylims(2));
% %         L = roundToNearest(max([x(:); y(:)]), 5, 'up');        
%         axis equal;
%         axis([0 L 0 L]); 
%     end    
end


%----------------------------------------------------------------------
function [tc_1, tc_2, idx] = rearrangeTrials(R, tc, ph_action, tr_action)
    [nOri, nSpf, nPh, nTrials] = size(R);     

    % 1. Reshape R into a matrix with each stimulus on a row, and trials on different columns     
    switch ph_action, 
        case 'average'
            R = mean(R,3);
            newdims = [nOri*nSpf,      nTrials];
        case 'max'
            R = max(R,[], 3);
            newdims = [nOri*nSpf,      nTrials];
        case 'itemize',
            newdims = [nOri*nSpf*nPh,  nTrials];
        case 'smooth',
            nFramesToAverage = deal(ph_action{2});
            nPh2 = ceil(nPh/nFramesToAverage);
            R2 = zeros(nOri, nSpf, nPh2, nTrials);
            for idx1 = 1:nPh2
                idx2 = (idx1-1)*nFramesToAverage + 1 : min(idx1*nFramesToAverage, nPh);
                R2(:,:,idx1,:) = mean(R(:,:,idx2,:), 3);
            end
%             [R, nPh] = deal(R2, nPh2);
            R = R2;
            newdims = [nOri*nSpf*nPh2,  nTrials];
    end
        
    R = reshape(R, newdims);
    
    
    % 2. use criteria to select some stimuli for the tuning curves
    select_method = tc{1};
    R_trial_av = nanmean(R,2);
    switch select_method        
        case {'ori', 'spf'}, [oriInds, spfInds] = deal(tc{[2,3]});  % only select oris (spfs) at a particular spf (ori)
%             nOriSpf = length(oriInds)*length(spfInds);
            idx = false(nOri, nSpf);
            idx(oriInds, spfInds) = true;
            idx = idx(:);
            
        case 'ori_spf',      [oriSpfInds] = tc{2};
%             nOriSpf = length(oriSpfInds);
            idx = oriSpfInds(:);
            
        case 'all_over',
            [op, x] = deal(tc{[2,3]});            
            if ~ibetween(x,0,1)
                error('threshold must be between 0 and 1');
            end
            idx = op(R_trial_av, x*max(R_trial_av));

        case 'all',
            [op, x] = deal(tc{[2,3]});            
            if ~ibetween(x,0,1)
                error('threshold must be between 0 and 1');
            end
            idx = op(R_trial_av, x*max(R_trial_av));
            
            
        case 'top'
            n = ceil(tc{2}*length(R_trial_av));
            idx = ord(R_trial_av, 'descend');
            idx = idx(1:n);
    end
    R = R(idx,:);
    
    
    % 3. Determine appropriate indices of trials, and divide data into 2 "tuning curves" 
    if nTrials > 1
        if (strcmp(tr_action, 'odd/even'))
            trialInds1 = 1:2:nTrials;
            trialInds2 = 2:2:nTrials;
        elseif (strcmp(tr_action, 'first/second'))
            trialInds1 = 1:nTrials/2;
            trialInds2 = nTrials/2:nTrials;
        end
    else
%         if strcmp(ph_action, 'itemize')
            warning('Not enough trials to compare across different trials: comparing across phases instead'); %#ok<WNTAG>
%             ph_action = 'compare';
%         end
        [trialInds1, trialInds2] = deal(1);
    end
    tc_1 = R(:,trialInds1);
    tc_2 = R(:,trialInds2);    

    % 4. Average over trials
%     tc_1a = mean(tc_1, 2);
%     tc_2a = mean(tc_2, 2);
    tc_1 = nanmean(tc_1, 2);
    tc_2 = nanmean(tc_2, 2);
    
end


