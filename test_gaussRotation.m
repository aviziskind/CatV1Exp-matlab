function test_gaussRotation

    N = 80;
    xs = linspace(-5, 5, N);
    ys = linspace(-5, 5, N);
    
    [xs_grid, ys_grid] = meshgrid(xs, ys);
    
    M = [4;1];    
    C = [1, .3; .3, 1];
    figure(1); clf;
    zs = gaussianN([xs_grid(:)'; ys_grid(:)'], M, C);
    zs = reshape(zs, N,N);
    
    M0 = [0;0];
    C0 = [.05, 0;0,.05];
    zC = gaussianN([xs_grid(:)'; ys_grid(:)'], M0, C0);
    zC = reshape(zC, N,N);
    
    zs = zs+zC/max(zC(:))*max(zs(:));
    
    
    h1 = surf(xs, ys, zs); 
    xlabel('x'); ylabel('y'); 
    hold on;
    view(2);
    axis equal square;
    3;
    h2 = surf(xs, ys, zs); 

    for th = linspace(0, 2*pi, 30)
        R = rotationMatrix(th);
        
        M2 = R*M;
        C2 = R*C*(R');
    
        zs2 = gaussianN([xs_grid(:)'; ys_grid(:)'], M2, C2);
        zs2 = reshape(zs2, N,N);

        set(h2, 'zdata', zs2);
        drawnow;
    end

end