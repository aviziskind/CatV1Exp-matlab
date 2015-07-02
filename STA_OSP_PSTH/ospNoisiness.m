function r = ospNoisiness(osp, w, omitFlag)

    global psthStatsSettings globW globOmitFlag

    collapseFunc = psthStatsSettings.ospPhCompressFcn;

    
    if ~exist('w', 'var') || isempty(w)
        if ~isempty(globW)
            w = globW;
        else
            w = .6;
        end
    end
    if ~exist('omitFlag', 'var') || isempty(omitFlag) && omitFlag
        if ~isempty(globOmitFlag)
            omitFlag = globOmitFlag;
        else
            omitFlag = false;
        end
        
    end
    
    osp = double(osp);
    if size(osp,3) > 1
        if strcmp(collapseFunc, 'max')
            osp = max(osp, [], 3);
        elseif strcmp(collapseFunc, 'mean')
            osp = mean(osp, 3);
        end
    end
    
    N = numel(osp);
    mn = min(osp(:));
    mx = max(osp(:));
    normOsp = (osp - mn)/ (mx - mn);
    
    % smooth OSP - in dim 1 (orientation), do circular smooth. in dim 2 (spf), do non-circular smooth. 
    smoothedOsp = gaussSmooth( gaussSmooth(normOsp, w, 1, 'circular', omitFlag), w, 2, [], omitFlag); %
    
%     r = sum( abs ( smoothedOsp(:) - normOsp(:) ) ./ (smoothedOsp(:) + .001) ) ;
    r = mean( abs ( smoothedOsp(:) - normOsp(:) )  );
    
end