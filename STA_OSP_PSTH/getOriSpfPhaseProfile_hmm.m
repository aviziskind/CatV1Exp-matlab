function [R, uori, usp, uph] = getOriSpfPhaseProfile_hmm(varargin)
    % input options: either 
    % (1)  Gid, cellId, spikeWindow and windowProfile (from which can calculate 'relContrOfFrameToSpike'
    %         ie. getOriSpfPhaseProfile(Gid, cellId, spikeWindow, windowProfile)
    %                   OR
    % (2) Gid (so know which stimulus type) and relContrOfFrameToSpike
    %         ie. getOriSpfPhaseProfile(Gid, relContrOfFrameToSpike)

    if (nargin == 2)
        [Gid, relContrOfFrameToSpike] = elements(varargin);        
    elseif (nargin == 4)
        [Gid, cellId, timeWindow, windowProfile] = elements(varargin);        
        relContrOfFrameToSpike = getParsedSpikes('frame',  Gid, cellId, timeWindow, windowProfile);
    else
        error('Syntax: call with (Gid, cellId, spikeWindow, windowProfile), or (Gid, relContrOfFrameToSpike)');
    end
    
    if iscell(relContrOfFrameToSpike)
        relContrOfFrameToSpike = [relContrOfFrameToSpike{:}];
    end

    [allOri_deg, allSp_pix, allPhase_deg, uori, usp, uph] = getOriSpPhaseForEachStimFrame(Gid);

    [frameStimIdsOS, oris, sps, phs] = getStimulusFrameSequence(Gid, 'OS');
        
    nOri = length(uori);
    nSp = length(usp);
    nPh  = length(uph);

    Theta = zeros(nOri, nSp, nPh);
    
    function y = h(osp_i, theta)
        [ori_i, sp_i, ph_i] = elements(osp_i);
        
        
        
    end
        

    function y = f(Q)
        
        
        y = exp ( )
    end
        
    
    R = zeros(nOri, nSp, nPh);
    nOccurencesOfStim = zeros(nOri, nSp, nPh);
    for iOri = 1:nOri
        curOris = allOri_deg == uori(iOri);
        for iSp = 1:nSp
            curOriSps = (curOris) & (allSp_pix == usp(iSp));
            for iPh = 1:nPh
                frmIdxs = (curOriSps) & (allPhase_deg == uph(iPh));
                nOccurencesOfStim(iOri, iSp, iPh) = length(frmIdxs);
                R(iOri, iSp, iPh) = sum( relContrOfFrameToSpike( frmIdxs ) );
            end
        end
    end
    % to do: rescale to # spikes / second (?)
    R = (R ./ nOccurencesOfStim) * median(nOccurencesOfStim(:));  % account for presentations that have extra frames.
    
end




%     oriIdxs = cell(1, nOri);
%     for iOri = 1:nOri
%         oriIdxs{iOri} = find(allOri_deg == uori(iOri));
%     end
%     spIdxs = cell(1, nSp);    
%     for iSp = 1:nSp
%         spIdxs{iSp} = find(allSp_pix == usp(iSp));
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
