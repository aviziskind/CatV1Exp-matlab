function [L_bin, R_bin, windowProfile] = getLRbins_windowprofile(timeWindow, PSTH)

    if isempty(PSTH)  % drifting gratings
        L_bin = 1;
        R_bin = 1;
        windowProfile = [];
        return
    end
    
    if strncmp(timeWindow, 'best', 4)
        if length(timeWindow) > 4
            offset_sign = switchh(timeWindow(5), {'M', 'P'}, [-1, 1]);
            offset = str2double(timeWindow(6))*offset_sign;
            timeWindow = 'best';
        else
            offset = 0;
        end
    end
        
    if ischar(timeWindow)
        switch timeWindow
            case 'best',  LR_bins = PSTH.timeWindow_bins + offset * [-1, 1];
            case 'stimw', LR_bins = PSTH.timeWindow_stimw;
        end
    elseif isnumeric(timeWindow)
        binCents = PSTH.bins;
        binEdges = binCent2edge(binCents);

        LR_bins(1) = indmin( abs(binEdges - timeWindow(1)) );
        LR_bins(2) = indmin( abs(binEdges - timeWindow(2)) ) - 1;  % this is the first edge of the next bin. want the previous bin.

%         %%
%         LR_bins(1) = find(PSTH.bins > timeWindow(1), 1);
%         LR_bins(2) = find(PSTH.bins > timeWindow(2), 1);
        
    end
    
    L_bin = LR_bins(1);
    R_bin = LR_bins(2);
    
    windowProfile = PSTH.vals(LR_bins(1) : LR_bins(2));
    3;

end