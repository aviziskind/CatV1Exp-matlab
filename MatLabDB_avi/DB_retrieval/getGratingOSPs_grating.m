function [ori_deg, sp_pix, phase_deg, uTempPeriods] = getGratingOSPs_grating(Gid)
        
    sd = siteDataFor(Gid);
    Did = sd.Did;
    presIds = sd.presIds;
    nPres = length(presIds);
%     nFramesPerPres = sd.stimulusInfo.nFramesPerPres;
    
    gratingOSP_path = [experimentDBPath 'Stimuli' filesep 'Grating' filesep];
    gratingOSP_name_mat = sprintf('GratingStimulus_Did_%d.mat', Did);

    filename = [gratingOSP_path gratingOSP_name_mat];
    if ~exist(filename, 'file')
        gratingParams = retrieveGratingPresParams(Did, presIds);        
        save(filename, '-struct', 'gratingParams');        
    else
        gratingParams = load(filename);
    end
    
    ori_deg = cell(1,nPres);
    sp_pix = cell(1,nPres);
    phase_deg = cell(1,nPres);
    tempPeriod_frm = zeros(1,nPres);
    
    for pres_i = 1:nPres        
        [ori_deg{pres_i}, sp_pix{pres_i}, phase_deg{pres_i}, tempPeriod_frm(pres_i)] = calcGratingOSPs_grating(gratingParams, pres_i);
    end
            
    [uTempPeriods, tp_lists] = uniqueList(tempPeriod_frm); 
    nTempPeriods = length(uTempPeriods);
            
    if nTempPeriods == 1
        ori_deg   = cat(1, ori_deg{:})';
        sp_pix    = cat(1, sp_pix{:})';
        phase_deg = cat(1, phase_deg{:})';
    elseif nTempPeriods > 1
        ori_deg   = cellfun(@(idx) cat(1, ori_deg{idx}),   tp_lists, 'un', 0);
        sp_pix    = cellfun(@(idx) cat(1, sp_pix{idx}),    tp_lists, 'un', 0);
        phase_deg = cellfun(@(idx) cat(1, phase_deg{idx}), tp_lists, 'un', 0);
    end
        
    %   ori_deg = mod(ori_deg, 360);
    %   uori = sort(unique(ori_deg));   % list of unique orientations
    %   usp = sort(unique(sp_pix));     % list of unique spatial frequencies
    %   uph = sort(unique(phase_deg));  % list of unique (spatial) phases

end

function [ori_deg, sp_pix, phase_deg, tempPeriod_frm] = calcGratingOSPs_grating(s, pres_i)

    [ori_deg, sp_pix, spatPhase_deg, tempPeriod_frm, tempPhase_deg, ...
        nDisplayedFrms, nPreBlankFrms, nPostBlankFrms, isStaticGrating] = deal(...
    s.ori_deg(pres_i), s.sp_pix(pres_i), s.spatPhase_deg(pres_i), s.tempPeriod_frm(pres_i), s.tempPhase_deg(pres_i), ...
        s.nDisplayedFrms(pres_i), s.nPreBlankFrms(pres_i), s.nPostBlankFrms(pres_i), s.isStaticGrating(pres_i));    

%     gratParams = struct('ori_deg', ori_deg, 'sp_pix', sp_pix, 'spatPhase_deg', spatPhase_deg, ...
%         'tempPeriod_frm', tempPeriod_frm, 'tempPhase_deg', tempPhase_deg, ...
%         'nDisplayedFrms', nDisplayedFrms, 'nPreBlankFrms', nPreBlankFrms, 'nPostBlankFrms', nPostBlankFrms, ...
%         'isStaticGrating', isStaticGrating);    

    nFramesInPres = nDisplayedFrms - (nPreBlankFrms + nPostBlankFrms);
    
    % ori_deg and sp_pix will all be the same value, whereas
    % phase_deg can either be a single value (for static (=flashed) gratings) or a 
    % different value for each frame (for drifting gratings), representing
    % the temporal phase of each frame.
        
    phase_deg_offset = mod(tempPhase_deg + spatPhase_deg, 360);
    if isStaticGrating  % static grating: phase is constant for entire presentation.
        phase_deg = phase_deg_offset;

    elseif ~isStaticGrating  % if drifting grating, phase is different for each frame                
        fracCycleOffset = phase_deg_offset/360;        
        fracsOfPres = (1./tempPeriod_frm)*[0:nFramesInPres-1];        
        fracsOfCycle = mod(fracsOfPres + fracCycleOffset, 1);
        
        if (tempPeriod_frm ~= round(tempPeriod_frm)) % for case where tempPeriod_frm is non-integer (= 43.737)
            phase_res = 1/tempPeriod_frm;
            fracsOfCycle = roundToNearest(fracsOfCycle, phase_res, 'down');
            fracsOfCycle = roundToDecimalPoint(fracsOfCycle, 5);
        end
        
        phase_deg = fracsOfCycle * 360;

        if max(abs(phase_deg-round(phase_deg))) < 1e-5  % for cases where have slight fractional differences in phase (eg Gid = 1048)
            phase_deg = round(phase_deg);
        end
        
        mn_phase = min(phase_deg);
        if mn_phase ~= 0;
            phase_deg = phase_deg-min(phase_deg); % for cases where nonzero phase_deg_offset would result in extra phases.  (eg Gid = 589)      
        end
    end    
    
    ori_deg = ori_deg(ones(size(phase_deg)));
    sp_pix  = sp_pix( ones(size(phase_deg)));    


end

function gratParams = retrieveGratingPresParams(Did, gratPresIds)
    hnd = dbOpenExpDb;
        
    dbParamFields = {'DBL_ORIENTATION_DEGR',  'DBL_SPATIAL_PERIOD_PIX', 'DBL_SPATIAL_PHASE_DEGR',   'DBL_TEMP_PERIOD_FRM',  'DBL_TEMP_PHASE_DEGR'};    
    dbFrameFields = {'LNG_N_DISPLAYED_FRM', 'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM', 'BLN_STATIC_GRATING'};

    [ori_deg, sp_pix, spatPhase_deg, tempPeriod_frm, tempPhase_deg, ...
        nDisplayedFrms, nPreBlankFrms, nPostBlankFrms, isStaticGrating, gratPresIds2] = ...        
    getFieldsFromDatabaseTable(hnd, [dbParamFields dbFrameFields, 'GRATING_PRES_ID'], 'TBL_GRATING_PRES', {'DATAFILE_ID', Did}, 'GRATING_PRES_ID');    
    assert(isequal(gratPresIds(:), gratPresIds2(:)));
    
    gratParams = struct('ori_deg', ori_deg, 'sp_pix', sp_pix, 'spatPhase_deg', spatPhase_deg, ...
            'tempPeriod_frm', tempPeriod_frm, 'tempPhase_deg', tempPhase_deg, ...
            'nDisplayedFrms', nDisplayedFrms, 'nPreBlankFrms', nPreBlankFrms, 'nPostBlankFrms', nPostBlankFrms, ...
            'isStaticGrating', isStaticGrating);    

end


%{
        tempPeriod_frm = getFieldsFromDatabaseTable(hnd, 'DBL_TEMP_PERIOD_FRM', 'TBL_GRATING_PRES', {'DATAFILE_ID', Did}, 'GRATING_PRES_ID');        
        [uTempPeriods, tp_lists] = uniqueList(tempPeriod_frm); nTempPeriods = length(uTempPeriods);
        [allOri_deg, allSp_pix, allPhase_deg] = deal(cell(1, nTempPeriods));
        
        if nTempPeriods > 1
            3;
        end
        for tp_i = 1:nTempPeriods         
            nFramesEachPres_tp_i = nFramesEachPres(tp_lists{tp_i});
            nTotalFrames_i = sum(nFramesEachPres_tp_i);
            cumSumFrames_i = [0; cumsum(nFramesEachPres_tp_i)];
            
            allOri_deg{tp_i}   = zeros(1, nTotalFrames_i);  
            allSp_pix{tp_i}    = zeros(1, nTotalFrames_i);
            allPhase_deg{tp_i} = zeros(1, nTotalFrames_i);
                        
            gratPresIds = getFieldsFromDatabaseTable(hnd, 'GRATING_PRES_ID', 'TBL_GRATING_PRES', {'DATAFILE_ID', Did; 'DBL_TEMP_PERIOD_FRM', uTempPeriods(tp_i)}, 'LNG_PRESENT_NO');
            for gi = 1:length(gratPresIds)
                [phase_deg, sp_pix, ori_deg] = getGratingOSPs_grating(hnd, gratPresIds(gi));
                allOri_deg{tp_i}(cumSumFrames_i(gi)+1:cumSumFrames_i(gi+1)) = ori_deg;
                allSp_pix{tp_i}(cumSumFrames_i(gi)+1:cumSumFrames_i(gi+1)) = sp_pix;
                allPhase_deg{tp_i}(cumSumFrames_i(gi)+1:cumSumFrames_i(gi+1)) = phase_deg;
            end
        end
        if nTempPeriods == 1
            allOri_deg   = allOri_deg{1};
            allSp_pix    = allSp_pix{1};
            allPhase_deg = allPhase_deg{1};            
        end
%}




% 
% plot(phase_deg(inds(1)+1:inds(2)), phase_deg2(inds(1)+1:inds(2)), 'b.', ...
%     phase_deg(inds(2)+1:inds(3)), phase_deg2(inds(2)+1:inds(3)), 'go', ...
%     phase_deg(inds(3)+1:inds(4)), phase_deg2(inds(3)+1:inds(4)), 'rs', ...
%     phase_deg(inds(4)+1:inds(5)), phase_deg2(inds(4)+1:inds(5)), 'm*', ...
%     phase_deg(inds(5)+1:inds(6)), phase_deg2(inds(5)+1:inds(6)), 'b.', ...
%     phase_deg(inds(6)+1:inds(7)), phase_deg2(inds(6)+1:inds(7)), 'go', ...
%     phase_deg(inds(7)+1:inds(8)), phase_deg2(inds(7)+1:inds(8)), 'rs', ...
%     phase_deg(inds(8)+1:inds(9)), phase_deg2(inds(8)+1:inds(9)), 'm*', ...
%     phase_deg(inds(9)+1:inds(10)), phase_deg2(inds(9)+1:inds(10)), 'b.', ...
%     phase_deg(inds(10)+1:inds(11)), phase_deg2(inds(10)+1:inds(11)), 'go')