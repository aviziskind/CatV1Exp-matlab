function [allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph, utp] = getOriSpPhaseForEachStimFrame(Gid)
            

    sd = siteDataFor(Gid);
    isFlashedGratingMovie = strncmp(sd.stimType, 'Movie:Flashed', 10);
    isDriftingGrating = strncmp(sd.stimType, 'Grating:Orientation', 12) || strncmp(sd.stimType, 'Grating:Spatial', 12);
        
    if isFlashedGratingMovie
        movieIds = sd.stimulusInfo.movieIds;
        
        [allOri_deg, allSp_pix, allPhase_deg] = getFlashGratingOSPs_fgmovie(movieIds);        
        uTempPeriods = 0;        
        
    elseif isDriftingGrating
        
        [allOri_deg, allSp_pix, allPhase_deg, uTempPeriods] = getGratingOSPs_grating(Gid);
        
    else
        error('Only call this function for Flashed / Drifting Grating stimuli, (not %s stimuli)', stimulusType);
    end

    
    if nargout > 3
        if ~iscell(allOri_deg)
            uori = unique(allOri_deg);
            usp =  unique(allSp_pix);
            uph  = unique(allPhase_deg);
            utp = S.uTempPeriods;
            
            if any(utp>0 && utp<50000) && (utp == round(utp))
                assert(length(uph) == utp);  % for drifting gratings -- make sure #of phases = temp_period in frame
            end
        else
            uori = cellfun(@unique, allOri_deg, 'un', 0);
            usp =  cellfun(@unique, allSp_pix, 'un', 0);
            uph  = cellfun(@unique, allPhase_deg, 'un', 0);
            utp =  uTempPeriods(:)';
        end
    end


end

%{
allGids = getAllGids;
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(allGids(i));
    progressBar(i);
end
getCorticalState('save');
%}


%     nFramesEachPres = getFieldsFromDatabaseTable(hnd, {'LNG_N_SUSTAINED_FRM'}, tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
%     if strcmp(stimulusType, 'Grating');
%         [nFadeIn, nFadeOut] = getFieldsFromDatabaseTable(hnd, {'LNG_N_FADEIN_FRM', 'LNG_N_FADEOUT_FRM'}, tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
%         nFramesEachPres = nFramesEachPres + nFadeIn + nFadeOut;
%     end


    %{
    hnd = dbOpenExpDb;
    Did = sd.Did;
        stimulusType = getStimulusTypeForDid(Did);
    tableName = getDatabaseTableForDid(Did, stimulusType);
    
    if ~exist('mode','var') || isempty(mode)
        mode = 'Shown';
    end
    switch mode
        case 'Planned'
            nSustainedFrms = getFieldsFromDatabaseTable(hnd, 'LNG_N_SUSTAINED_FRM', tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
            if dbDoesFieldExist(hnd, 'LNG_N_FADEIN_FRM', tableName)
                [nFadeInFrms, nFadeOutFrms] = getFieldsFromDatabaseTable(hnd, {'LNG_N_FADEIN_FRM', 'LNG_N_FADEOUT_FRM'}, tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
            else
                [nFadeInFrms, nFadeOutFrms] = deal(0);
            end
            nFramesEachPres = nSustainedFrms + (nFadeInFrms + nFadeOutFrms);
        case 'Shown'
            [nDisplayedFrms, nPreBlankFrms, nPostBlankFrms] = getFieldsFromDatabaseTable(hnd, {'LNG_N_DISPLAYED_FRM', 'LNG_N_PRE_BLANK_FRM', 'LNG_N_POST_BLANK_FRM'}, tableName, {'DATAFILE_ID', Did}, 'LNG_PRESENT_NO');
            nFramesEachPres = nDisplayedFrms - (nPreBlankFrms + nPostBlankFrms);    
    end
    %}

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


%{
allDgGids = getAllGids('d'); n = length(allDgGids);
progressBar('init-', n)
for i = 1:n
    getOriSpPhaseForEachStimFrame(allDgGids(i)) 
    progressBar(i);
end


%}


