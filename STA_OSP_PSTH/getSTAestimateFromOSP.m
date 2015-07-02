function STA = getSTAestimateFromOSP(Gid, OSP)

    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid);
    Did = dbLookup('Did',  'Gid', Gid);
    hnd = dbOpenExpDb;
    [ncols, nrows] = getFieldsFromDatabaseTable(hnd, {'LNG_N_COLUMNS', 'LNG_N_ROWS'}, 'TBL_MOVIE_PRES', {'DATAFILE_ID', Did}, [], 1, 1);
    
    oris = sort(unique(allOri_deg));
    sps = sort(unique(allSp_pix));
    phs  = sort(unique(allPhase_deg));

    %take top  5% of responses
    [bestOriSpPhs, bestOriSpPhsInds] = findBestOriSpPhs(OSP, oris, sps, phs);
    frameDims = [ncols, nrows];
    STA = zeros(frameDims);
    getMovieStimulusFrame('load',  'Gid', Gid);
    for i = 1:size(bestOriSpPhs,1)
%         [ori_deg, sp_pix, ph_deg] = elements( bestOriSpPhs(i,:) );
        weight = sub2indV( size(OSP), bestOriSpPhsInds(i,:) );
        stimulatoryFrameInds = findFrameIndsWithCertainOSP(Gid, bestOriSpPhs(i,:));
        frm = getMovieStimulusFrame(stimulatoryFrameInds(1));

        STA = STA + weight*frm;
        %generateGratingFrame(frameDims, deg2rad(ori_deg), sp_pix, ph_deg );
    end


end
