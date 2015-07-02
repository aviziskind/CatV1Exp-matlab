function [ori, sp] = findHighestCommonOriSps(OSPs, n)

    nOSPs = length(OSPs);
    for oi = 1:nOSPs
        thisOSP = OSPs{oi};
        [osp, oris, sps, phs] = elements(thisOSP);
        subplot(1,nOSPs, oi);
        imagesc3({oris, sps, phs}, osp, 3, {'Orientation', 'Sp period', 'Sp Phase'});
        colorbar;
    end


end