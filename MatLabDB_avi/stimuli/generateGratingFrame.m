function frame = generateGratingFrame(frameDims, ori_deg, spatPeriod_pix, spatPhase_deg,  tempPeriod_sec, tempPhase_deg, t_sec, isSquare )

    % dimensions: [k] = 1/pix;  [omega] = 1/t, [theta] = radians, 
    %             [phi] = radians, [t] = sec.
    
    if ~exist('spatPhase_deg', 'var') || isempty(spatPhase_deg)
        spatPhase_deg = 0;
    end
    if ~exist('tempPeriod_sec', 'var') || isempty(tempPeriod_sec)
        tempPeriod_sec = 0;
    end
    if ~exist('tempPhase_deg', 'var') || isempty(tempPhase_deg)
        tempPhase_deg = 0;
    end    
    if ~exist('t_sec', 'var') || isempty(t_sec)
        t_sec = 0;
    end
    isSquare = exist('isSquare', 'var') && ~isempty(isSquare) && isSquare;
    
    isStatic = (tempPeriod_sec == 0);
    if ~(isStatic)  % if drifting grating, ignore spatPhase parameter
        spatPhase_deg = 0;
    end
    
    if isSquare
        waveFunction = @(phi) sign(sin(phi)); % square wave.
    else
        waveFunction = @sin;
    end
    
    if exist('frameDims', 'var') && ~isempty(frameDims);
        Nx = frameDims(1);
        Ny = frameDims(2);        
    else
        [Nx, Ny] = deal(1024, 768);  % full field grating
    end
    xs = 0:Nx-1;
    ys = 0:Ny-1;
    
    theta = deg2rad(ori_deg);
    k   = 2*pi/spatPeriod_pix;
    phi_x = deg2rad(spatPhase_deg);

    omega = 1/tempPeriod_sec;
    phi_t = deg2rad(tempPhase_deg);
    
    [xs_grid, ys_grid] = meshgrid(xs, ys);
    rotatedXY = rotationMatrix( theta ) * [xs_grid(:), ys_grid(:)]';
    xs_grid_rot = reshape(rotatedXY(1,:), size(xs_grid));    

    if isStatic
        frame = waveFunction(k * xs_grid_rot + phi_x + phi_t)';
    else
        frame = waveFunction(k * xs_grid_rot + phi_x  - omega*t_sec  + phi_t)';
    end
    
end

%     zs  = Amin + (Amax-Amin) * waveFunction(k * xs_grid - omega*t + phi);    
%     figure(10); imagesc(zs);        colormap('gray'); colorbar; axis equal square tight xy; 
%     figure(11); imagesc(zsRotated); colormap('gray'); colorbar; axis equal square tight xy; 




%     function y = squareWave2(phi)
%           % same as squareWave(phi) (this is a quicker implementation in C, since
%           % it skips the sine function), but it's a little slower in matlab 
%         y = zeros(size(phi));
%         dphi = mod(phi, 2*pi);
%         up = ((dphi < 0) & (dphi < -pi)) | ((dphi > 0) & (dphi < pi));
%         down = ~up;
%         y(up) = 1;
%         y(down) = -1;
%     end

%     frame = generateGratingFrame(ori_deg, spatPeriod_pix, spatPhase_deg, tempPeriod_sec, tempPhase_deg,  isSquare, t_sec, frameDims)