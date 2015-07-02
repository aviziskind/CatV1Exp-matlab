function [dcR, corR, oris, sps, phs] = getOriSpfPhaseProfile_decor(Gid, varargin)
    % input options: either 
    % (1)  Gid, spkTsRelToFrame_ms, timeWindow_ms, frameLength_ms)  (much faster) 
    % (1)  Gid, cellId (and do everything from scratch)

    if (nargin == 2)
        [Gid, cellId] = elements(varargin);    
        
        [PSTHdata, spkTsRelToFrame_ms] = getPSTHforCell(Gid, cellId, isFlashGratingStimulus); warning('dc:getPSTHformat', 'no longer correct format');
        frameLength_ms = PSTHdata.frameLength_ms;
        timeWindow_ms = PSTHdata.timeWindow;                
        
    elseif (nargin > 2)
        [spkTsRelToFrame_ms, timeWindow_ms] = elements(varargin);
        if nargin == 4            
            frameLength_ms = varargin{3};                
        else
            frameLength_ms = getFrameLength('Gid', Gid);
        end
    end
    
    [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms);
    nBinsPerExtFrame = nBinsPerFrame * nFramesPerExtFrame;
    
    % Generate PSTHs of every stimulus (*not* grouping by spatial phase)
    [frameStimIdsOSP, oris, sps, phs] = getStimulusFrameSequence(Gid, 'OSP');
      uStimIdsOSP = unique(frameStimIdsOSP);
    PSTHs = zeros(nBinsPerExtFrame,  length(uStimIdsOSP));
    for stim_i = uStimIdsOSP(:)'
        frameIndsForStimI = find(frameStimIdsOSP == stim_i);
        [PSTHbins, PSTHs(:,stim_i)] = calcPSTH( spkTsRelToFrame_ms(frameIndsForStimI), extFrameLength_ms, length(frameIndsForStimI), [], nBinsPerExtFrame);        
    end
    
    % 2. Find which frames are included in the time window;
    windowBinInds(1) = find( PSTHbins > timeWindow_ms(1), 1, 'first');
    windowBinInds(2) = find( PSTHbins < timeWindow_ms(2), 1, 'last' );
    windowFramesInds = ceil(windowBinInds / nBinsPerFrame);
    windowFrames = windowFramesInds(1):windowFramesInds(2);
    
    % 5. Decorrelating the PSTHs.
    dcPSTHs = decorrelatePSTHs(PSTHs, frameStimIdsOSP, nBinsPerFrame, nFramesPerExtFrame, windowFrames);
        
    nOri = length(oris);
    nSp = length(sps);
    nPh  = length(phs);

    corR = reshape(sum(PSTHs, 1),   [nOri, nSp, nPh]);    
    dcR  = reshape(sum(dcPSTHs, 1), [nOri, nSp, nPh]);    
    
end




%  tic;
%     R2 = zeros(nOri, nSp, nPh);
%     for iOri = 1:nOri
%         for iSp = 1:nSp
%             for iPh = 1:nPh
%                 frmIdxs = (allOri_deg == oris(iOri)) & (allSp_pix == sps(iSp)) & (allPhase_deg == phs(iPh));
%                 R2(iOri, iSp, iPh) = sum( relContrOfFrameToSpike( frmIdxs ) );
%             end
%         end
%     end
%     toc;