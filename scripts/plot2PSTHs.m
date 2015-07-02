function plot2PSTHs(binC, vals1, vals2)

    doSmooth = false;
    binE = binCent2edge(binC);
    binE_plot = binE(1:end-1);
        
    y1 = vals1(:); 
    if doSmooth        
%         y1 = extt(y1);
        binE_itp = linspace(binE(1), binE(end), 100);
        y1_itp = interp1(binE, y1, binE_itp, 'spline');        
        plot(binE_itp, y1_itp);
    else
        stairs(binE_plot, y1, 'linewidth', 2);             
    end    
    
    if nargin > 2
        hold on;
        y2 = vals2(:);
        if doSmooth           
%             y2 = extt( y2 );
            y2_itp = interp1(binE, y2, binE_itp, 'spline');        
            plot(binE_itp, y2_itp, 'r');
            axis tight;
        else
            wh_pix = getObjDims(gca, 'pixels');
            x_range = binE(end)-binE(1);
            d = x_range / wh_pix(1);
        	stairs(binE_plot+d*3, y2, 'r', 'linewidth', 2);
        end
    end
        
    
    xlim([binE(1), binE(end)]);

end