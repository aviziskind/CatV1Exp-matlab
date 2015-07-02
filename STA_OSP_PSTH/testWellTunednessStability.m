function testWellTunednessStability

    S = load('flashGratingData8.mat');
    allOSPs = S.allOSPs;
    
    osp_i = 9;
    [groupId, cellId, oris, sps, phs, OSP] = elements( allOSPs(osp_i) );

    
    function OSPnew = switchAround(OSP)
%        find top 10
        n = 30;
        [nx, ny] = size(OSP);
        [tmp, bestInds] = maxElements(OSP, n);
        OSPnew = OSP;
        
        % try randomly switching
        for i = 1:n
            [x,y] = elements( bestInds(i,:) );
            [rx, ry] = elements(  1-floor(rand(2,1)*3)  ); % [-1, 0, 1];
            if (x == 1) && (rx == -1) || (x == nx) && (rx == 1)
                rx = -rx;
            end
            if (y == 1) && (ry == -1) || (y == ny) && (ry == 1)
                ry = -ry;
            end
            oldInd = sub2ind([nx ny], x,y);
            newInd = sub2ind([nx ny], x+rx, y+ry);
            OSPnew([oldInd, newInd]) = OSPnew([newInd, oldInd]);
        end
    end
        
    
%     dw = 100;
    OSP = sum(OSP, 3);

    w = findHowWellTunedOSP(OSP);
      
    while true
        j = 1;
        OSPnew = switchAround(OSP);
        assert( abs( sum(OSPnew(:)) - sum(OSP(:))) < 1e-5 );
        w_new = findHowWellTunedOSP(OSPnew);
            
        if w_new > w;
            fprintf('%d) successful swap (%f vs %f) \n', j, w, w_new );
            OSP = OSPnew;
            w = w_new;
            imagesc(OSP);
            title(['w = ' num2str(w)]);
            drawnow;
        else
            fprintf('%d) no swap (%f vs %f) \n', j, w, w_new );
        end
        
        j = j+1;            
    end

        
        
        
end
    