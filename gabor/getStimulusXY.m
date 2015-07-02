function [xs_deg, ys_deg] = getStimulusXY(si, dsampFactor)
    if isnumeric(si)
        Gid = si;
        sd = siteDataFor(Gid);
        si = sd.stimulusInfo;
    end

    [origin_x_pix, origin_y_pix, pixPerBlock, nrows, ncols, screenH_pix, screenW_pix, degPerBlock] = deal(...
    si.origin_x, si.origin_y, si.pixPerBlock, si.nrows, si.ncols, si.screenH_pix, si.screenW_pix, si.degreesPerBlock);
   
    if (nargin < 2) || isempty(dsampFactor)
        dsampFactor = 1;
    end

    xs_pix_UL = origin_x_pix + [0:dsampFactor:ncols-1]*pixPerBlock;
    ys_pix_UL = origin_y_pix + [0:dsampFactor:nrows-1]*pixPerBlock;

    xs_pix = xs_pix_UL - screenW_pix/2;
    ys_pix = ys_pix_UL - screenH_pix/2;
    
    degPerPix = degPerBlock/pixPerBlock;
    xs_deg = xs_pix*degPerPix;
    ys_deg = ys_pix*degPerPix;
    
    show = 0;
    if show        
        center_xy_pix = [screenW_pix/2;  screenH_pix/2];            
        rect_pix = [xs_pix_UL(1), ys_pix_UL(1), xs_pix(end)-xs_pix(1), ys_pix(end)-ys_pix(1)];        
        ax_pix = [1 screenW_pix, 1 screenH_pix];

        figure(456); clf; hold on; box on;
        axis ij;
        axis(ax_pix);        
        rectangle('position', rect_pix);        
        plot(center_xy_pix(1), center_xy_pix(2), 'r*')
        title('pixels (from upper-left corner of screen)');
        
        ax_deg = [ax_pix(1:2) - center_xy_pix(1), ax_pix(3:4) - center_xy_pix(2)]*degPerPix;
        rect_deg = [xs_deg(1), ys_deg(1), xs_deg(end)-xs_deg(1), ys_deg(end)-ys_deg(1)];
        
        figure(457); clf; hold on; box on;
        axis ij;        
        axis(ax_deg);        
        rectangle('position', rect_deg);
        plot(0, 0, 'r*');
        title('degrees (from center of screen)');
    end
    3;
        
        
end
