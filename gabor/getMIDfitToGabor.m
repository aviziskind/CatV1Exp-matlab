function [MID_fit, xs, ys] = getMIDfitToGabor(arg1, gaborParams)
    if length(arg1) == 1
        Gid = arg1;
        [xs, ys] = getStimulusXY(Gid);
        [x_grid, y_grid] = meshgrid(xs, ys);
        XX = [x_grid(:), y_grid(:)];
        rf_size = [length(xs), length(ys)];
    else
        XX = arg1;
        L = sqrt(length(XX));
        rf_size = L*[1, 1];
        [xs, ys] = deal(1:L);
    end

    MID_fit = reshape(gaborP(gaborParams, XX), rf_size);        
end