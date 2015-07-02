function STA = getSTAfromOSP(Gid, OSP, scaling)
    if nargin < 3
        scaling = 'aroundZero';
    end

    frameStimIds = getStimulusFrameSequence(Gid);    
    if size(OSP, 4) > 1 % if multiple trials.
        OSP = mean(OSP, 4);
    end    
    nStim = numel(OSP);
    
    [uStimIds, idx_firstOccurence] = unique(frameStimIds, 'first');
    stim_firstFrameIds = reshape(idx_firstOccurence, size(OSP));    
    
    getFrame = getFrameRetrieverFunction(Gid);    
    getFrame('load', Gid, 'scaling', scaling);
    frameDims = getFrame('size');
    STA = zeros(frameDims);
    
    % rescale OSP
    OSP = double(OSP);
%     OSP = OSP/mean(OSP(:));  % rescaling makes even/odd-trial STAs not sum to all-trial STAs
    
    for stim_i = 1:nStim
        wgt = OSP(stim_i);
        frm_id = stim_firstFrameIds(stim_i);         
        frm = getFrame(frm_id);

        STA = STA + (wgt*frm)/nStim;
        %generateGratingFrame(frameDims, deg2rad(ori_deg), sp_pix, ph_deg );
    end
    getFrame('close');

end
