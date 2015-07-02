function tmp_plotStdMID(Gid, cellId)
[mid_fileName] = mid_getPreferredMIDfile(Gid, cellId);
S = load(mid_fileName);
MID = S.MID;
v_MID = S.v_MID{1};
nx = size(MID,1);
MID_std = std( v_MID, [], 3);

figure(55); clf;
for i = 1:4
    h_ax(i) = subplot(2,4,i); imagesc(v_MID(:,:,i)); axis square; 
end

h_ax(5) = subplot(2,4,5:6);
imagesc(MID); axis square; colorbar;
3;
matchAxes('C', h_ax(1:5))
3;

v_MID_vec = reshape(v_MID, [nx*nx, 4]);
c_mtx = corr(v_MID_vec);
[~, c_vec] = pearsonRm(v_MID_vec);

title(sprintf('jcc = %.2f', mean(c_vec)));

% show = 'dist';
% show = 'std';
% show = 'MID_corr';
show = 'none';
    
switch show
    case 'MID_corr'

        nPix = 3;
        [corr_mtx, corr_X] = pixelCorr(S.v_MID{1}, nPix);

        h_ax(6) = subplot(2,4,7:8);
        imagesc(corr_mtx); axis square
        caxis([-1 1]);

        imagesc(corr_X); axis square
        caxis([-1 1]);
        colorbar;
        
        title('Dot Product Matrix');
        axes(h_ax(5)); hold on;
        contour(corr_X);
        
    case 'std';
        h_ax(7) = subplot(2,4,7:8);
%         a = mean(abs(MID(:)))/2;
%         MID_std_scl = MID_std ./ (a + abs(MID));
        imagesc( MID_std )
        axis square
        colorbar
        title('Std Dev');
        
    case 'dist',
        h_ax(7) = subplot(2,4,7:8);
        dist_mtx = pixelOffDiag(v_MID);
        imagesc( dist_mtx );
        axis square
        
        
end

%}

        axes(h_ax(1));  %#ok<MAXES>
        title(sprintf('[%d, %d]', Gid, cellId));


3;

3;
end


function [corr_mtx, corr_sizeX] = pixelCorr(X, npix)
    [M,N,nJack] = size(X);
    assert(M == N);    
    
    wind_size = 2*npix+1;
    nc = N-wind_size+1;
    dpix = [-npix:npix];
    corr_mtx = zeros(nc, nc);    
    corr_sizeX = zeros(N);    
    V = zeros(wind_size^2, nJack);
    
    idx_lower = find ( tril(ones(nJack), -1) );
    
%     pix_corr = zeros(nJack);
    
    i0_idxs = npix+1:M-npix;
    j0_idxs = npix+1:N-npix;
    
    for i = i0_idxs
        for j = j0_idxs
            for jack_i = 1:nJack
                v_i = X(i + dpix, j + dpix, jack_i);
                V(:,jack_i) = v_i(:);
            end
            
            
            pix_corr = zeros(nJack);
            for ci = 1:nJack
                for cj = 1:ci
%                     pix_corr(ci,cj) = pearsonR(V(:,ci), V(:,cj));                                        
                    pix_corr(ci,cj) = normDotProd(V(:,ci), V(:,cj));
                end
            end
%             [pix_corr, pix_corr_vec] = pearsonRm(V);
            3;
            
%             pix_corr2 = corr(V);            
%             corr_sizeX(i,j) = mean( pix_corr_vec );             
            corr_sizeX(i,j) = mean( pix_corr(idx_lower) );              %#ok<FNDSB>
            
        end
    end

    
    3;
    if nargout > 1
%         corr_sizeX = zeros(N,N);
        for i = 1:N
            i_idx = bound(i-npix, 1, nc)+npix;
            
            for j = 1:N
                j_idx = bound(j-npix, 1, nc)+npix;
                
                if ((i <= npix) || (i >= N-npix+1)) || ...
                        ((j <= npix) || (j >= N-npix+1))
                    assert(corr_sizeX(i,j) == 0)
                    assert(corr_sizeX(i_idx,j_idx) ~= 0)
                    
                    corr_sizeX(i,j) = corr_sizeX(i_idx, j_idx);
                    
                end
                
%                 corr_sizeX(i,j) = corr_mtx(i_idx, j_idx);
            end
        end                
    end      
    
    corr_mtx = corr_sizeX(i0_idxs, j0_idxs);
    3;
    

end

function dist_mtx = pixelOffDiag(X)
    [M,N,nJack] = size(X);
    assert(M == N);    
    
    dist_mtx = zeros(N);
    
    H = 1/4*ones(4);
    for i = 1:N
        for j = 1:N
            x = X(i,j,:); x = x(:);
            xp = H*x;
            
            d = norm(x-xp);
            dist_mtx(i,j) = d;            
        end        
    end
    
end