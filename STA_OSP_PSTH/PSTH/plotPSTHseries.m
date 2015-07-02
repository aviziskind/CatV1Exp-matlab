function hs = plotPSTHseries(PSTHseries, plotStyle, plotOptions, plotLabels)
    % plotOptions is a cell with 2 optional values: {stimOrder, [plotParams]}
    % for 2D plots, plotParams = [figStart, nrows, ncols];
    % for 3D plots, plotParams = [figureId];

    nSeries = length(PSTHseries);
    nStims = size(PSTHseries{1},2);
    
    % default for stimOrder
    if (nargin < 3) || isempty(plotOptions{1})
        stimOrder = 1:nStims;
    else
        stimOrder = plotOptions{1};
    end
    
    % default for plotParams
    if (nargin < 3) || length(plotOptions) == 1 || isempty(plotOptions{2})
        if strcmp(plotStyle, '2D')
            nrowsDefault = iff(nSeries > 2, nSeries, nSeries*2);
            plotParams = [520, nrowsDefault, 5];         
        elseif strcmp(plotStyle, '3D')
            plotParams = [520];
        end
    else
        plotParams = plotOptions{2};
    end
    
    
    if strcmp(plotStyle, '2D')
        [figStart, nrows, ncols] = elements(plotParams);
        maxFigs = 20;
        gridSubPlot(nrows,ncols, [figStart maxFigs], nSeries);
        hs = zeros(nSeries, nStims);
        for ser_i = 1:nSeries
            PSTHs = PSTHseries{ser_i};
            if isempty(PSTHs), break; end;

            for stim_i = 1:nStims                    
                lastOne = gridSubPlot(ser_i); 
                hs(ser_i, stim_i) = bar(PSTHs(:,stimOrder(stim_i)), 1, color(ser_i));
                set(gca, 'Xtick', []);                    
                if (ser_i == 1)
                    title(num2str(stim_i)); 
                end
                if lastOne, break; end;
            end
            
%             matchAxes('Y', nonzeros(hs(ser_i,:)) );
        end
        matchAxes('Y', nonzeros(hs([1,3],:)) );
    elseif strcmp(plotStyle, '3D')
        figId = plotParams(1);
        
        if ~isempty(figId) && (figId > 0)
            figure(figId); 
        else
            figure;
        end
        hs = zeros(1,nSeries);
        nrows = floor(sqrt(nSeries));
        ncols = ceil(nSeries / nrows);                    
        for ser_i = 1:nSeries
            PSTHs = PSTHseries{ser_i};
            subplot(nrows, ncols, ser_i);
            if isempty(PSTHs), break; end;
            
            h = bar3(PSTHs(:,stimOrder));
            hs(ser_i) = get(h(1), 'Parent');
            set(gca, 'ydir', 'normal'); view(94, 5);
            axis tight;
            if exist('plotLabels', 'var')
                title(plotLabels{ser_i});
            end
            
        end
        
    end
            
end