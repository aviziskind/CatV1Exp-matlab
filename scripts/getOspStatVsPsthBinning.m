function [stats, r, r_odd, r_even] = getOspStatVsPsthBinning(Gid, psthBins, allPsthVals, l_bin, r_bin_orig, shiftForEvenNBins, curPsth, outputVars, CellId)

%     global tmpVals;
    shiftForEvenNBins = isempty(shiftForEvenNBins) || shiftForEvenNBins; % assume don't shift unless specified.    
    storeOsps = exist('CellId', 'var') && ~isempty(CellId);
    
    global psthStatsSettings
    
%     nudge0pval_amt = psthStatsSettings.nudge0pval_amt;        
    pvalfunc = iff(psthStatsSettings.storeLogOfPvalues, @(x) negLogPval_loc(x, psthStatsSettings.nudge0pval_amt), ...
                                                        @(x) x);    
    setClass = str2func(psthStatsSettings.statsPrecision);
                                                    
%     [rep_pval, rep_tstat, rep_offset, r_odd, r_even] 
%     okOutputs = {'rep_p', 'rep_tstat', 'rep_slope', 'rep_offset',  'rho', 'rho_p', 'rho_t',   'tau', 'tau_p', 'tau_t', 'r_var', 'r_entropy'};
    okOutputs = {'cc', 'cc_p', 'cc_t',  ...
                'rho', 'rho_p', 'rho_t',  'rho_nz', 'rho_p_nz', 'rho_t_nz',  'rho_p_nznz', 'rho_t_nznz', ...
                'tau', 'tau_p', 'tau_t', ...
                'rep_p', 'rep_tstat', 'rep_slope', 'rep_offset',  ...
                'r_var', 'r_entropy', 'noisiness', 'n_unique'};
            
    for i = 1:length(outputVars)
        if ~any(strcmp(outputVars{i}, okOutputs))
            error('Invalid outputVar: %s', outputVars{i})
        end
    end
        
%     strCcmp
    doPearsonR = any( strncmp('cc', outputVars, 2) );
    doKendallTau = any( strncmp('tau', outputVars, 3) );
    doSpearmanRho = any( strcmp('rho', outputVars) ) || any( strcmp('rho_p', outputVars) );
    doSpearmanRho_nz = any( strcmp('rho_nz', outputVars) ) || any( strcmp('rho_p_nz', outputVars) );
    doSpearmanRho_nznz = any( strcmp('rho_nznz', outputVars) ) || any( strcmp('rho_p_nznz', outputVars) );

    doRegressionTest = any( strncmp('rep', outputVars, 3) );
    doVariance = any( strcmp('r_var', outputVars) );
    doEntropy = any( strcmp('r_entropy', outputVars) );
    doNoisiness = any( strcmp('noisiness', outputVars) );
    doNUnique   = any( strcmp('n_unique', outputVars) );
    
%     getEvenOddTrials = doKendallTau || doSpearmanRho || doSpearmanRho_nz || doSpearmanRho_nznz || doRegressionTest || nargout >= 2;
%     getSingleOSP = doVariance || doEntropy || doNoisiness || doNUnique || nargout >= 4;
        
%     [dims, repType, frmLength_ms] = parseStimType(stimType);
    
    

    nbins_orig = r_bin_orig-l_bin+1;
    shiftCenter = shiftForEvenNBins && ~odd(nbins_orig);
    
    if shiftCenter && (r_bin_orig == size(allPsthVals,1))
        stats = nan(1,length(outputVars));
        return;
    end
    
    if storeOsps        
        [r, r_oe] = getOspDataForPsthWindow(Gid, CellId, psthBins, allPsthVals, l_bin, r_bin_orig, curPsth, {'osp', 'osp_oe'}); % this function calls 'calc', but also saves the results
    else
        [r, r_oe] = calcOspForPsthWindow(Gid, {psthBins, allPsthVals}, l_bin, r_bin_orig, shiftForEvenNBins, curPsth, {'osp', 'osp_oe'});
    end
    r_odd = r_oe(:,:,:,1);
    r_even = r_oe(:,:,:,2);
    r = r(:); r_odd = r_odd(:); r_even = r_even(:); 

    
    function [stat, stat_p, stat_t] = doCorr(x, y, varargin)
         if isempty(x) || isempty(y) || ~any(x) || ~any(y)
             [stat, stat_p, stat_t] = deal(nan);
         else
             try 
                 [stat, stat_p, stat_t] = myCorr(x, y, varargin{:});
%                  [stat2, stat_p2] = corr(x, y, varargin{:});
%                  assert(stat == stat2)
%                  assert(stat_p == stat_p2);
             catch ME
                 keyboard;
                 beep;
                 [stat, stat_p, stat_t] = deal(nan);
             end
         end
         stat_p = pvalfunc(stat_p);
         
         stat   = setClass(stat);
         stat_p = setClass(stat_p);
         stat_t = setClass(stat_t);
    end

    
    
    if doPearsonR
        [cc, cc_p, cc_t] = doCorr(r_odd, r_even, 'type', 'pearson', 'tail', 'gt');      %#ok<NASGU>
%         [cc, cc_p] = myCorr(r_odd, r_even, 'type', 'pearson', 'tail', 'gt');    
    end

    
    if doKendallTau
%         [tau, tau_p] = corr(r_odd, r_even, 'type', 'kendall', 'tail', 'gt');        
        [tau, tau_p, tau_t] = doCorr(r_odd, r_even, 'type', 'kendall', 'tail', 'gt');      %#ok<NASGU>    
%         [tau2, tau_p2] = myCorr(r_odd, r_even, 'type', 'kendall', 'tail', 'gt');
%         assert(tau == tau2);
%         assert(tau_p == tau_p2);
    end

    
    if doSpearmanRho
        [rho, rho_p, rho_t] = doCorr(r_odd, r_even, 'type', 'spearman', 'tail', 'gt');        %#ok<NASGU>  
%         [rho2, rho_p2] = myCorr(r_odd, r_even, 'type', 'spearman', 'tail', 'gt');        
%         assert(rho== rho2);
%         assert(rho_p == rho_p2);
    end

    if doSpearmanRho_nz
        idx_eitherNonZero = (r_odd | r_even);
        r_odd_nz  = r_odd(idx_eitherNonZero);
        r_even_nz = r_even(idx_eitherNonZero);        
        [rho_nz, rho_p_nz, rho_t_nz] = doCorr(r_odd_nz, r_even_nz, 'type', 'spearman', 'tail', 'gt');   %#ok<NASGU>      
    end

    if doSpearmanRho_nznz
        idx_bothNonZero = (r_odd & r_even);
        r_odd_nznz  = r_odd(idx_bothNonZero);
        r_even_nznz = r_even(idx_bothNonZero);        
        [rho_nznz, rho_p_nznz, rho_t_nznz] = doCorr(r_odd_nznz, r_even_nznz, 'type', 'spearman', 'tail', 'gt'); %#ok<NASGU>
    end
    
    if doRegressionTest
        show = [];
        getTstatFlag = iff(any(strcmp('rep_tstat', outputVars)), 1, []);
        [rep_tmp, rep_p, rep_slope, rep_offset] = regressionSlopeTtest(r_odd, r_even, .01, '+', show, getTstatFlag);     %#ok<NASGU>
        rep_p = pvalfunc(rep_p);       
        if getTstatFlag
            rep_tstat = rep_slope;     
        end
        rep_p = setClass(rep_p);          %#ok<NASGU>
        rep_slope = setClass(rep_slope);  %#ok<NASGU>
        rep_tstat = setClass(rep_tstat);  %#ok<NASGU>
    end

    
    if doVariance
        r_var =          -var(r / sum(r));             %#ok<NASGU>
    end    

    if doEntropy
        r_entropy = sum ( entropy (r / sum(r)));     %#ok<NASGU>
    end    
        
    if doNoisiness        
        noisiness = ospNoisiness( r );     %#ok<NASGU>
    end
    
    if doNUnique
        n_unique = length(unique( r(:) ));     %#ok<NASGU>
    end
    
    stats = cellfun(@eval, outputVars);
                
    
end



function p = negLogPval_loc(p, amt)
    if p == 0
        p = -log10(amt);
    else
        p = -log10(p);
    end
    % if p was 1, it will now be 0. want to keep 0 for 'unassigned',
    % though, so assign 0's to -0.01
    if p == 0
        p = cast(-0.1, class(p));
    end      
    
end
%         x = r_odd(:); y = r_even(:);    
%         X = [ones(size(x)), x];
%         b = regress(y, X);        
%         rep_pval = b(2); % first output argument is thus actually the slope.
    
    
%         binOSP_odd  = zeros([dims, nBins]);
%         binOSP_even = zeros([dims, nBins]);
%         nStim = prod(dims);
%         for bin_i = 1:nBins
%             binOffset = dims * (bin_i-1);
%             binOSP_odd( [1:nStim] + binOffset) = psths_odd(bin_i,:) ;
%             binOSP_even( [1:nStim] + binOffset) = psths_even(bin_i,:) ;
%         end            
%         


%{
%     shiftSecondWindowIf2x8 = false;    

    shiftAmount_ms = 100 + frmLength_ms/2;

%     shiftSecondWindow = shiftSecondWindowIf2x8; %&& is2x8;
%     if shiftSecondWindow
%         zeroBin_idx = find(psthBins >= 0, 1);
%         if ~(r_bin >= zeroBin_idx);            
%             shiftAmount_bins = ms2bin(shiftAmount_ms);
%             psth_val_idx2 = psth_val_idx2 + shiftAmount_bins;
%             if any(psth_val_idx2 >= zeroBin_idx)
%                 psthWindowStart = psthBins(1)-binW/2;
%                 psth_val_idx2 = psth_val_idx2 + ms2bin(shiftAmount_ms+psthWindowStart);
%             end
%         end        
%         
%     end


%}

