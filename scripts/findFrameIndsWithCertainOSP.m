function stimulatoryFrameInds = findFrameIndsWithCertainOSP(Gid, oriSpfPhs)

    [allOri_deg, allSp_pix, allPhase_deg] = getOriSpPhaseForEachStimFrame(Gid);
    OSP_matrix = [allOri_deg(:), allSp_pix(:), allPhase_deg(:)];    
    n = size(oriSpfPhs,1);
    allFrameInds = cell(1,n);
    for i = 1:n
        allFrameInds{i} = findRows(oriSpfPhs(i,:), OSP_matrix);
    end
    stimulatoryFrameInds = sort([allFrameInds{:}]);

end