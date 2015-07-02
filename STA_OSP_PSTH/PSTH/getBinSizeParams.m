function [nBinsPerFrame, nFramesPerExtFrame, extFrameLength_ms] = getBinSizeParams(frameLength_ms, extFrameLength_ms)

    global idealBinSize_ms idealExtFrameLength_ms;

    if isempty(idealBinSize_ms)
        idealBinSize_ms = 4+(1/6); %
    end
    if isempty(idealExtFrameLength_ms)
        idealExtFrameLength_ms = 150;
    end
    
    if nargin < 2
        extFrameLength_ms = idealExtFrameLength_ms;        
    end    
    nFramesPerExtFrame = ceil(extFrameLength_ms/frameLength_ms);
    extFrameLength_ms = nFramesPerExtFrame * frameLength_ms;

    nBinsPerFrame = round(frameLength_ms / idealBinSize_ms);            
end