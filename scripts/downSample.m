function zs_coarse = downSample(zs, n)

%     figure(23);
    [nx, ny, nSamp] = size(zs);
    nx_new = ceil(nx/n);
    ny_new = ceil(ny/n);
    zs_coarse = zeros(nx_new, ny_new, nSamp);

    for i = 1:nx_new
        old_inds_i = (i-1)*n+1:i*n;
        if old_inds_i(end) > nx
            old_inds_i(old_inds_i > nx) = [];
        end
        for j = 1:ny_new
            old_inds_j = (j-1)*n+1:j*n;
            if old_inds_j(end) > ny
                old_inds_j(old_inds_j > ny) = [];
            end
            zs_select = zs( old_inds_i, old_inds_j, : );    
            if nSamp > 1
                zs_coarse(i,j, :) = mean(mean( zs_select, 1 ),2);
            else
                zs_coarse(i,j) = mean(zs_select(:));
            end
        end
    end

end