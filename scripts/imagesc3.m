function [h_ax, h_im] = imagesc3(ax, data, collapseDim, collapseMode, axlabels)

    if ~exist('collapseDim', 'var') || isempty(collapseDim)
        collapseDim = 3;
    end
    
    collapsedDim = [];
    bSingleDims = (size(data) == 1);    
    if nnz(bSingleDims) == 0  % only collapse if there is not already a singleton dimension
        collapsedDim = collapseDim;
        
        switch collapseMode
            case 'mean', data = nanmean(data, collapseDim);        
            case 'max',  data = max(data, [], collapseDim);        
            case 'pref',
                [tmp, inds_max] = maxElement(data);                     % ie.   data = data(inds(1), :, :)
                idx = arrayfun(@(n) 1:n, size(data), 'Un', false);      %   or  data = data(:, inds(2), :)
                idx{collapseDim} = inds_max(collapseDim);               %   or  data = data(:, :, inds(3))
%                 fprintf('pick out row # %d of dimension %d\n', inds_max(collapseDim), collapseDim)
                data = data(idx{:});
%             otherwise
%                 idx = arrayfun(@(n) 1:n, size(data), 'Un', false);      %   or  data = data(:, inds(2), :)
%                 idx{collapseDim} = inds_max(collapseDim);
                
                
        end
        
    end    
    data = squeeze(data);   
    singleDims = find(bSingleDims);    
    nonSingleDims = find(~bSingleDims);
    
    if ~isempty(collapsedDim)
        nonSingleDims = setdiff(nonSingleDims, collapsedDim);
    end
    
    if isvector(data)  % only 1 non-singleton dimension
        
        h_im = imagesc(ax{nonSingleDims(1)}, ax{singleDims(1)}, data');        
        if exist('axlabels', 'var') && ~isempty(axlabels)
            xlabel(axlabels(nonSingleDims(1)));
        end
                
    else 
        
        axInd1 = nonSingleDims(1);
        if length(nonSingleDims) > 1
            axInd2 = nonSingleDims(2);
        else
            axInd2 = singleDims(1);
        end
        h_im = imagesc(ax{axInd2}, ax{axInd1}, data);

        if exist('axlabels', 'var') && ~isempty(axlabels)
            xlabel(axlabels{axInd2});
            ylabel(axlabels{axInd1});
        end
    end
    h_ax = gca;

    [ny, nx] = size(data);
    xblockSize = diff(xlim)/nx;
    yblockSize = diff(ylim)/ny;        
%     if xblockSize >
    if ny > 1.5*nx
        xy_ratio = .7;
    elseif nx > 1.5*ny
        xy_ratio = 1/(.7);
    else
        xy_ratio = 1;
    end;
        
    daspect([xy_ratio * xblockSize, yblockSize, 1]);
    
end

