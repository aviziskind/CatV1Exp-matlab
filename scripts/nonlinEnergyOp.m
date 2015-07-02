function x = nonlinEnergyOp(x, dim, bartlett_N)
    sizeX = size(x);
    
    if length(sizeX) > 2
        error('Not implemented for more than 2 dimensions');
    end

    if nargin < 2
        dim = find(sizeX > 1, 1);
    end

    idx_mid = 2:(size(x, dim)-1);
    switch dim
        case 1,   x(idx_mid,:) = x(idx_mid,:).^2 - x(idx_mid-1,:).*x(idx_mid+1,:);
        case 2,   x(:,idx_mid) = x(:,idx_mid).^2 - x(:,idx_mid-1).*x(:,idx_mid+1);
    end        
    
    if (nargin >= 3) && (bartlett_N > 0)
        bartlett_window = bartlett(bartlett_N);
        bartlett_window = bartlett_window/sum(bartlett_window);
                
        switch dim
            case 1,
                for j = 1:sizeX(2)
                    x(:,j) = conv(x(:,j), bartlett_window, 'same');
                end
                
            case 2,
                for j = 1:sizeX(1)
                    x(j,:) = conv(x(j,:), bartlett_window, 'same');
                end
                
        end
    end       
        
        

end


%         case 1,   x(idx_mid,:,:) = x(idx_mid,:,:).^2 - x(idx_mid-1,:,:).*x(idx_mid+1,:,:);
%         case 2,   x(:,idx_mid,:) = x(:,idx_mid,:).^2 - x(:,idx_mid-1,:).*x(:,idx_mid+1,:);
%         case 3,   x(:,:,idx_mid) = x(:,:,idx_mid).^2 - x(:,:,idx_mid-1).*x(:,:,idx_mid+1);
