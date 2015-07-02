function OSP = predictOSPfromSTA(Gid, STA)
    
%     P_s_v
%     P_s_v_spk
    

    [frameStimIds, uOri, uSpf, uPh] = getStimulusFrameSequence(Gid);            
    
    [uStimIds, idx_firstOccurence] = unique(frameStimIds, 'first');    
    nStim = length(uStimIds);
    stim_firstFrameIds = idx_firstOccurence; %reshape(idx_firstOccurence);    
    
    getFrame = getFrameRetrieverFunction(Gid);    
    getFrame('load', Gid, 'scaling', 'aroundZero');    
    
    OSP = zeros(length(uOri), length(uSpf), length(uPh));
    STA = double(STA);
    
    for stim_i = 1:nStim
        frm_id = stim_firstFrameIds(stim_i);         
        frm = double( getFrame(frm_id) );
        
        OSP(stim_i) = sum(STA(:) .* frm(:));
        %generateGratingFrame(frameDims, deg2rad(ori_deg), sp_pix, ph_deg );
    end    
    getFrame('close');



end