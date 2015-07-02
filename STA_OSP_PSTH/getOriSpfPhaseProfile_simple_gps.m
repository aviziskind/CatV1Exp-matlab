function [R, R_full, f, uori, usp, uph, utp] = getOriSpfPhaseProfile_simple_gps(varargin)

    limitToPres = [];    
    % input options: either 
    % (1)  Gid, cellId, spikeWindow and windowProfile (from which can calculate 'relContrOfFrameToSpike'
    %         ie. getOriSpfPhaseProfile(Gid, cellId, spikeWindow, windowProfile, [, meanFiringRate]] )
    %                   OR
    % (2) Gid (so know which stimulus type) and relContrOfFrameToSpike
    %         ie. getOriSpfPhaseProfile(Gid, relContrOfFrameToSpike, [ [bckgRate, meanFiringRate] ])

    switch nargin
        case {2, 3}
            [Gid, relContrOfFrameToSpike] = elements(varargin);
            if nargin == 3 && ~isempty(varargin{3});
                rates = varargin{3};
            end
            
        case {4, 5}
            [Gid, cellId, timeWindow, windowProfile] = elements(varargin);
            [relContrOfFrameToSpike] = getParsedSpikes('frame',  Gid, cellId, timeWindow, windowProfile);
            if nargin == 5 && ~isempty(varargin{5});
                rates = varargin{5};
            end
        otherwise            
            error('Syntax: call with (Gid, cellId, spikeWindow, windowProfile, [bckg, [av]]), or (Gid, relContrOfFrameToSpike, [bckg, [av]])');
    end
   
    [uOri, uSp, uPh, nPermsOrCycles, nReps] = dbGetUniqueOriSpPh('Gid', Gid);
    nRepsTotal = nPermsOrCycles * nReps;
    if nRepsTotal == 1
        error('stim:singleTrial', 'Only 1 repetition of each grating stimulus was presented (need at least 2)');
    end

    if exist('rates', 'var')
%         bckgRate = rates{1};
%         if length(rates) > 1
%             meanFiringRate = rates{2};
%         end
    end
    
    if iscell(relContrOfFrameToSpike)
        relContrOfFrameToSpike = [relContrOfFrameToSpike{:}];
    end
    
    function [R, R_full, f, uori, usp, uph] = getOSP_profile(allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph)
    
        nOri = length(uori);
        nSp = length(usp);
        nPh  = length(uph);
        nStimTypes = nOri*nSp*nPh;
        nTrialsAv = round(length(allOri_deg)/nStimTypes);

        if ~isempty(limitToPres)
            hnd = dbOpenExpDb; Did = dbLookup('Did',  'Gid', Gid);
            stimTableName = getDatabaseTableForDid(Did);
            nSustFrames = getFieldsFromDatabaseTable(hnd, 'LNG_N_SUSTAINED_FRM', stimTableName, {'DATAFILE_ID', Did});
            cumFrames = [0; cumsum(nSustFrames)];

            framesInPres = cumFrames(limitToPres)+1:cumFrames(limitToPres+1);
            framesInOtherPres = setdiff(1:cumFrames(end), framesInPres);
        end

        % don't know what the max # of trials is across all stimuli 
        R_full = nan(nOri, nSp, nPh, nTrialsAv);  
%         R      = zeros(nOri, nSp, nPh, 'single');
        nTrialsTheseParams = zeros(nOri, nSp, nPh); 
        nTrialsMax = nTrialsAv;
        
        for iOri = 1:nOri
            curOris = allOri_deg == uori(iOri);
            for iSp = 1:nSp
                curOriSps = (curOris) & (allSp_pix == usp(iSp));
                for iPh = 1:nPh
                    frmIdxs = (curOriSps) & (allPhase_deg == uph(iPh));                
                    if ~isempty(limitToPres)
                        frmIdxs(framesInOtherPres) = 0;
                    end
                    n = nnz(frmIdxs);
                    if n > nTrialsMax
                        R_full(:,:,:,nTrialsMax+1:n) = NaN;  % expand with NaNs, not with zeros.
                        nTrialsMax = n;
                    end                    
                    nTrialsTheseParams(iOri, iSp, iPh) = n;
                    R_full(iOri, iSp, iPh, 1:n) = relContrOfFrameToSpike( frmIdxs );
%                     if n > 0
%                         R(iOri, iSp, iPh) = nanmean(R_full(iOri, iSp, iPh, 1:n));
%                     end
                end
            end
        end
        
        if (flashedOrDrifting(Gid) == 2)  % for drifting gratings - have to discard first cycle which may contain artifacts (of sudden change from no stimulus)
            nCyclesToSkip = 1;
            nCycles = ceil(nPermsOrCycles);
            okTrialsInds = bsxfun(@plus, [0:nReps-1]'*nCycles, 1+nCyclesToSkip:nCycles)';
            okTrialsInds = okTrialsInds(:)';
        else
            okTrialsInds = 1:nTrialsMax;
        end
        assert(nTrialsMax == max(nTrialsTheseParams(:)));    
        
        R = nanmean(R_full(:,:,:,okTrialsInds), 4);
        R(nTrialsTheseParams == 0) = 0;

        
%         assert(max(abs(R(:) - R2(:))) < .01);
%         R = R2;
        
        % rescale to # spikes / second (?)
        R_full = R_full / frameLength_sec; % convert from # eff spikes --> #spikes/sec.
        R = R           / frameLength_sec; % convert from # eff spikes --> #spikes/sec.
        
        
        f.scale = 255/max(R_full(:));
        f.idxnan = uint32( find(isnan(R_full(:))) );
        R_full = uint8(R_full*f.scale);

%         figure; imagesc(nanmean(nanmean(R_full, 3),4));
%         fprintf('%.3f\n', max(R_full(:)));
%         R_full_orig = R_full;


        % for consistency with "recalculations" later on, do the lossy compression once now:        
%         R_full = single(R_full)/f.scale;
%         R_full(f.idxnan) = nan;
%         R = nanmean(R_full(:,:,:,okTrialsInds), 4);        
        
    end    
        
    frameLength_sec = getFrameLength('Gid', Gid, 'sec');
    [allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph, utp] = getOriSpPhaseForEachStimFrame(Gid);
    
    if length(utp) == 1
        [R, R_full, f, uori, usp, uph] = getOSP_profile(allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph);
    else
        [R, R_full, f, uori, usp, uph] = cellfun(@getOSP_profile, allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph, 'un', 0);
    end
    
%     if exist('firingRate', 'var')
%         meanRate2 = mean(R(:));
%         discrep = (meanRate2 - meanFiringRate)/meanFiringRate * 100;
%         fprintf('OSP: Check:  mean rate 1: %2.3f. mean rate 2: %.3f  (discrepancy: of %.4f %%)\n', meanFiringRate, meanRate2, discrep);
%     end
%     R = (R ./ nOccurencesOfEachParamSet) * median(nOccurencesOfEachParamSet(:));  % account for presentations that have extra frames.           
        
%     figure(132);
%     ss = struct('R', single(R), 'ori', uori, 'sp', usp, 'ph', uph);
%     imageOSP(ss);
%     3;
end




%     oriIdxs = cell(1, nOri);
%     for iOri = 1:nOri
%         oriIdxs{iOri} = find(allOri_deg == uori(iOri));
%     end
%     spIdxs = cell(1, nSp);    
%     for iSp = 1:nSp
%         spIdxs{iSp} = find(allspf_pix == usp(iSp));
%     end
%     phIdxs = cell(1, nPh);    
%     for iPh = 1:nPh
%         phIdxs{iPh} = find(allPhase_deg == uph(iPh));
%     end
%     
%     R = zeros(nOri, nSp, nPh);
%     for iOri = 1:nOri
%         for iSp = 1:nSp
%             for iPh = 1:nPh
%                 frmIdxs = intersectq(oriIdxs{iOri}, spIdxs{iSp}, phIdxs{iPh});
%                 R(iOri, iSp, iPh) = sum( relContrOfFrameToSpike( frmIdxs ) );
%             end
%         end
%     end



%     for p = 1:nPh
%         subplot(2, nPh, p );       imagesc(mean(R_full(:,:,p,trialInds1), 4));
%         set(gca, 'xtick', [], 'ytick', []);
%         subplot(2, nPh, p + nPh ); imagesc(mean(R_full(:,:,p,trialInds2), 4));
%         set(gca, 'xtick', [], 'ytick', []);
%     end
% 
%    



%         if nSp > 1            
%             % get spatial frequency
%             spfTuningCurve = mean( R(ori_peak_ind, :,:) ,3); % at peak ori, sum across phases
%             nullSpfCurve = ones(size(spfTuningCurve)) * sum(spfTuningCurve)/length(spfTuningCurve);
%             
%             if (length(nonnans(spfTuningCurve)) > 1)
%                                 
%                 [spf_sel,spf_pval] = histChiSqrTest(spfTuningCurve, nullSpfCurve, alpha);
%                 if showWorking
%                     figure(103); clf;
%                     bar(spfTuningCurve); hold on;
%                     xlim(.5+[0, length(spfTuningCurve)])
%                     drawHorizontalLine( nullSpfCurve(1), 'color', 'r');
%                     title( sprintf('Spf non-uniform: h = %d. p = %2.4g', spf_sel,spf_pval));
%                 end
% %                 uniformCDF = [min(spfTuningCurve), max(spfTuningCurve); 0 1]';
% %                 [spf_sel, spf_pval] = kstest( spfTuningCurve, uniformCDF, alpha );            
%             end
%     %         [spf_sel_h2, spf_pval2] = lillietest( spfTuningCurve-mean(spfTuningCurve), alpha );
%         end
%         stats.spatialfrequencySelective = spf_sel; 
%         stats.spatialfrequencySelectivePval = spf_pval;    


%             spfTuningCurves_MS = bsxfun(@minus, spfTuningCurves, mean(spfTuningCurves,1));
%                 figure(105); clf;
%                 hist(rs_shuffled_MS);
%                 drawVerticalLine(r_actual_MS, 'color', 'r');
%                 drawVerticalLine(m2, 'color', 'g');
%                 drawVerticalLine([m2 + s2,m2 - s2], 'color', 'g','linestyle', ':');
%                 drawVerticalLine(median(rs_shuffled_MS), 'color', 'b');
%                 title( sprintf('p = %2.4g', p2));                            
