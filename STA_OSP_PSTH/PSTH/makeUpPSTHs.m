function r = makeUpPSTHs(nBinsPerExtFrame, nStimuli, binWindow, nStd)
    if nargin < 4, 
        nStd = 1;
    end
    rMax = 15;
    
    % make up some r values    
    r = nStd * randn(nBinsPerExtFrame, nStimuli) * 1 ; %noise

    m = mean(binWindow);
    s = diff(binWindow)/2;

    binInds = binWindow(1):binWindow(2);
    %add some peaks
    for i = 1:nStimuli
        thisPeak = histc( s*randn(1,5000)+m, 1:nBinsPerExtFrame )';
        thisPeak = thisPeak / max(thisPeak) * rMax * ((nStimuli-i)/nStimuli)^3;
        r(binInds,i) = r(binInds,i) + thisPeak(binInds);
    end
    r = r(:, ord(sum(r,1), 'descend') );   % scramble the order

end
