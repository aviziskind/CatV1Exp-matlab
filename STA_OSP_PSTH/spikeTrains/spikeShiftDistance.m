function dist = spikeShiftDistance(spk1, spk2, T, spk1_weights, fig_id, ax)

    dbug = exist('fig_id', 'var');
    
    if exist('T', 'var') && ~isempty(T)
        % allow for periodic boundary conditions:
        spk2 = sort(mod(spk2, T));
        spk2  = [spk2-T; spk2; spk2+T];
    end


    inds = binarySearch(spk2, spk1, 1, 2);
    dists =  abs( spk1 - spk2(inds) ).^2 ;
    if exist('spk1_weights', 'var') && ~isempty(spk1_weights)
        dists  = dists .* spk1_weights;
    end
    dist = sum(dists);    
    
    
    if dbug
        figure(fig_id); clf; hold on;
        X1 = [1; 1] * spk1(:)';     Y1 = (1 + [-0.5; 0.5]) * ones(1,length(spk1));
        if exist('spk1_weights', 'var') && ~isempty(spk1_weights)
            Y1(1,:) = 1 - [0.5 * spk1_weights];
        end
        
        X2 = [1; 1] * spk2(:)';      Y2 = (2 + [-0.5; 0.5]) * ones(1,length(spk2));
        line( X1, Y1, 'Color', 'g' );  % spk1
        line( X2, Y2, 'Color', 'r' );  % spk2
        quiver(X1(1,:), mean(Y1,1), X2(1,inds)-X1(1,:), mean(Y2(:,inds),1)-mean(Y1,1), 'AutoScale', 'off');
        if exist('ax', 'var')
            xlim(ax);
        else
            xlim( [ min([spk1(:); spk2(:)]), max([spk1(:); spk2(:)]) ] );
        end
        axis ij;
        ylim([0.5 2.5]);
        hold off;
    end
        

end

   