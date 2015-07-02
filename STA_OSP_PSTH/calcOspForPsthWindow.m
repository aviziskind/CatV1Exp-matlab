function varargout = calcOspForPsthWindow(Gid, arg2, l_bin, r_bin_orig, shiftForEvenNBins, curPsth, whatToCalculate, getHistArgs, meanRate)
    % default: [r, r_odd, r_even] = ...
    %    if whatToCalculate = 'phase':
    % [r_ph, r_odd_ph, r_even_ph] = ...
    
    global psthStatsSettings
    compressMethod = psthStatsSettings.ospPhCompressFcn; % options: 'mean'/'max';        
    ospPhCompress = switchh(compressMethod, {'max', 'mean'}, {@(osp) max(osp, [], 3), @(osp) nanmean(osp, 3) } );    

    % get bins & psth values
    if isnumeric(arg2) && (length(arg2) == 1)         
        cellId = arg2;
        if ~exist('getHistArgs', 'var') || isempty(getHistArgs)
            getHistArgs = {};
        end
        
        [psthBins, allPsthVals, meanRate] = dbGetCellSpkStimHists(Gid, cellId, struct('trialGrouping', 'individual', getHistArgs{:})); % get with all trials separately
    elseif iscell(arg2)
        psthBins    = arg2{1};
        allPsthVals = arg2{2};        
    end        
    [nOri, nSp, nPh, nTrials, isCph] = getGratingStimType(Gid);
    
    if ~exist('meanRate', 'var')
        meanRate = mean(allPsthVals(:));
    end
    
    rescaleToMatchMeanRate = 1;% && exist('meanRate', 'var');
    
    canCalculate = {'osp',      'osp_oe',  ...
                    'osp_ph',   'osp_ph_oe', ...
                    'osp_full', 'meanRate'};
    
%     if ~exist('whatToCalculate', 'var') || isempty(whatToCalculate)
%         whatToCalculate = {'osp', 'osp_odd', 'osp_even'};
%     elseif ischar(whatToCalculate) && strcmp(whatToCalculate, 'phase')
%         whatToCalculate = {'osp_ph', 'osp_ph_odd', 'osp_ph_even'};
%     end
    if ischar(whatToCalculate)
        whatToCalculate = {whatToCalculate};
    end

    if ~all( strCcmp(whatToCalculate, canCalculate))
        error('Provided unknown variables to calculate');
    end
    if nargout ~= length(whatToCalculate)
        warning('Mismatch of # outputs with number of variables to calculate');
    end
    
    returnWithoutPhase = any(strCcmp(whatToCalculate, {'osp', 'osp_oe'}));
    returnWithPhase = any(strCcmp(whatToCalculate, {'osp_ph', 'osp_ph_oe'}));
    returnWithAllTrials = any(strCcmp(whatToCalculate, 'osp_full'));
        
    if returnWithAllTrials && (~isCph && (size(allPsthVals, 3) <= 2))
        error('Can''t return osp_full if don''t have all trials');
    end
    
    % determine precise window size
    if isempty(shiftForEvenNBins)
        shiftForEvenNBins = 1;
    end
    
    nbins_orig = r_bin_orig-l_bin+1;
    shiftCenter = shiftForEvenNBins && ~odd(nbins_orig);
    
    nbins = nbins_orig;    
    r_bin = r_bin_orig;
    
    if shiftCenter && (r_bin < size(allPsthVals,1))
        nbins = nbins+1;
        shift_wgt = ones(nbins, 1);
        shift_wgt([1,end]) = .5;
        r_bin = r_bin+1;
    end
    psth_val_idx = l_bin:r_bin;
    if length(psth_val_idx) == length(curPsth)  % if provided just the relevant part of the psth (within the window)
        curPsth_val_idx = 1:length(curPsth);
    elseif length(psth_val_idx) == length(psthBins) % if provided the entire psth (and want us to now take the relevant part)
        curPsth_val_idx = psth_val_idx;
    end
    
    haveCurPsth = exist('curPsth', 'var') && ~isempty(curPsth);    
    
    
    % weights (for psth or for center shift.)
    wgts = ones(nbins, 1);    
    if haveCurPsth
        curPsth = curPsth(:);
        wgts = wgts .* curPsth(curPsth_val_idx(:));
    end
    if shiftCenter
        wgts = wgts .* shift_wgt;
    end          
%     wgts = wgts(:) / sum(wgts);

    
    
    % allPsthVals_oe should be of size [# of bins in PSTH x nTotalStimuli x 2] 
    % (odd trials in first page, even trials in second page).    
    nTrials = size(allPsthVals, 3);
    haveAllTrials = nTrials > 2;
    if haveAllTrials        
        idx_odd = 1:2:nTrials;
        idx_even = 2:2:nTrials;
        oe_factor = length(idx_even)/length(idx_odd);
        allPsthVals_oe = cat(3, mean(allPsthVals(:,:,idx_odd),  3), ...
                                mean(allPsthVals(:,:,idx_even), 3)*oe_factor );
    else
        allPsthVals_oe = allPsthVals;        
    end                                
    
    
    % each of these three cell arrays has 3: 
    %  (1) mean (all trials), (2) odd trials, (3) even trials 
    psths_oe = cell(1,3);
    osps     = cell(1,3);
    osps_ph  = cell(1,3);
        
    psths_oe{1} =    mean(allPsthVals_oe(psth_val_idx, :, :), 3);   % 1 - average over all trials     
    psths_oe{2} =         allPsthVals_oe(psth_val_idx, :, 1);       % 2 - odd trials
    psths_oe{3} =         allPsthVals_oe(psth_val_idx, :, 2);       % 3 - even trials
    
    [rescaleFactor_osp_ph, rescaleFactor_osp] = deal(1);
    
    for i = 1:length(psths_oe)
        if (haveCurPsth || shiftCenter);
            psths_oe{i}  = bsxfun(@times, psths_oe{i},  wgts(:));
        end
        allR  = nanmean(psths_oe{i},  1);        
        osp_ph_i = reshape(allR, [nOri, nSp, nPh]);
        
        if returnWithPhase  % return full 36x10x8 profile                        
            if i == 1 && rescaleToMatchMeanRate  % use same rescaling for odd & even (based on all-trial scaling)
                rescaleFactor_osp_ph = meanRate / mean(osp_ph_i(:)); 
            end
            
            osps_ph{i} = osp_ph_i * rescaleFactor_osp_ph;
        end
        
        if returnWithoutPhase % return 36x10 profile (compress along phase dimension)
            if strcmp(compressMethod, 'max')
                if (i == 1)            
                    osp_i    = ospPhCompress(osp_ph_i);        
                    select_idx = getSelectMaxIdx(osp_ph_i);
                else            
                    osp_i = osp_ph_i(select_idx);
                end
            elseif strcmp(compressMethod, 'mean')
                osp_i  = ospPhCompress(osp_ph_i);           
            end
            
            if i == 1 && rescaleToMatchMeanRate
                rescaleFactor_osp = meanRate / mean(osp_i(:)); 
            end
            
            osps{i} = osp_i * rescaleFactor_osp;            
        end
        
    end

    
    if returnWithAllTrials
        psths_full =  allPsthVals(psth_val_idx, :, :);   
        psths_full =  bsxfun(@times, psths_full,  wgts(:));
        allR  = nanmean(psths_full,  1);        
        osp_full = reshape(allR, [nOri, nSp, nPh, nTrials]);  
        
        osp_full = osp_full / mean(osp_full(:))*meanRate; %#ok<NASGU>  % scale to match mean firing rate.
    end    
    
        
    if returnWithoutPhase
        osp      = osps{1};         
        assert(~any(isnan(osp(:))));
        if any(strcmp(whatToCalculate, 'osp_oe'))
            osp_odd  = osps{2};
            osp_even = osps{3};
            osp_oe   = cat(4, osp_odd, osp_even); %#ok<NASGU>
        end
    end
    
    if returnWithPhase
        osp_ph      = osps_ph{1};   %#ok<NASGU>
        if any(strcmp(whatToCalculate, 'osp_ph_oe'))
            osp_ph_odd  = osps_ph{2};   
            osp_ph_even = osps_ph{3};   
            osp_ph_oe   = cat(4, osp_ph_odd, osp_ph_even); %#ok<NASGU>
        end
    end    
    
    for i = 1:length(whatToCalculate)
        varargout{i} = eval(whatToCalculate{i});
    end
            

end



function select_idx = getSelectMaxIdx(osp_ph)
    useF1phaseForAmbiguousCases = false;

    [max_ph, phmax_idx] = max(osp_ph, [], 3);
    
    if useF1phaseForAmbiguousCases    
        [dbl_max_ori, dbl_max_sp] = find ( (max_ph > 0) & sum(bsxfun(@eq, max_ph, osp_ph), 3) > 1 ); 
        nPh = size(osp_ph, 3);
        phases = linspace(0, 360, nPh+1); phases = phases(1:nPh);    
        for i = 1:length(dbl_max_ori)
            oi = dbl_max_ori(i);
            si = dbl_max_sp(i);
            phase_tc = osp_ph(oi, si, :);

            ph_max_rad = getF1phase(phases, phase_tc(:) ) ;
            if isnan(ph_max_rad), continue,  end

            ph_max_idx = round( ph_max_rad / (2*pi)*nPh )+1;
            if ph_max_idx == nPh+1, ph_max_idx = 1; end

            phmax_idx(oi, si) = ph_max_idx;
            if isnan(phmax_idx(oi, si))
                3;
            end
        end
    end
    
    ori_sp_offset = reshape(1:360, 36, 10);
    ph_offset     = (phmax_idx-1)*360;
    select_idx = ori_sp_offset + ph_offset;
    
end


% 
%             okTrialsInds = bsxfun(@plus, [0:nRep-1]'*nCycles, 2:nCycles)'; 
%             okTrialsInds = okTrialsInds(:)';
%         
%             oddTrials_idx = odd(1:length(okTrialsInds));
%         
%             R_odd  = nanmean(R_full(:,:,okTrialsInds( oddTrials_idx)), 3);
%             R_even = nanmean(R_full(:,:,okTrialsInds(~oddTrials_idx)), 3);



%     splitPSTHwindowIfCph = true;
%         cphSplitTime_ms = 70;
%     ...
%     cphWindow = (size(allPsthVals,3) == 4);
%     if cphWindow && splitPSTHwindowIfCph
%         bins_before =  ( psthBins <  cphSplitTime_ms );
%         bins_after  =  ( psthBins >= cphSplitTime_ms );                
%         allPsthVals = cat(1, allPsthVals(bins_before, :, [1 2], :), ...
%                              allPsthVals(bins_after,  :, [3 4], :) );
%     end                  
