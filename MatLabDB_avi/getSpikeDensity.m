function [ts,fs] = getSpikeDensity(spikeTimes, windowSize)

    plotResults = true;

    halfWindowSize = windowSize/2;    
    N = 10; % N points per window

    t_first = spikeTimes(1);
    t_last  = spikeTimes(end);

    dt = windowSize / N;
    ts = t_first-dt : dt : t_last+dt;
    nPoints = length(ts);
    fs = zeros(1,nPoints);    
    
    
    Ws = [max( ts(:)-halfWindowSize, t_first),  min( ts(:)+halfWindowSize, t_last)];

    for i = 1:nPoints
        fs(i) = length(elementsInRange(spikeTimes, Ws)) / diff(W);
    end
    

%     for i = 1:nPoints
%         progressBar(i);
%         W = [max(ts(i)-halfWindowSize, t_first),  min(ts(i)+halfWindowSize, t_last)];
%         fs(i) = length(elementsInRange(spikeTimes, W)) / diff(W);
%     end

    if plotResults
        plot(ts,fs);
    end

end