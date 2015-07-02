function cycleThroughFigures(hs, numTimes, delay_sec, titles)
    % syntax:
    %    alternateBetweenFigures(h1, h2, ... hn);
    %    alternateBetweenFigures(..., nTimes, delay);

    numFigs = length(hs);
    if (nargin < 2)
        numTimes = 3;
    end
    if (nargin < 3)
        delay_sec = 0.25;
    end
    
    for fj = 1:numFigs
        set(hs(fj), 'Visible', 'on');
    end
    axis auto;
    v = axis;
    
    for ti = 1:numTimes

        for fi = 1:numFigs
            set(hs(fi),  'Visible', 'on');
            for fj = setdiff(1:numFigs, fi)
                set(hs(fj), 'Visible', 'off');
            end
            axis(v);
            title(titles{fi});
            drawnow;
            pause(delay_sec);
        end
        
    end
end
