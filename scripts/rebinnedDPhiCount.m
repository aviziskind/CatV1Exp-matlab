    function cnt_binned = rebinnedDPhiCount(dPhis, dPhiNull_Count, newBinEdges)
%         if gratingType == 1
%             cnt_binned  = dPhiNull_Count;
%             return;
%         end        
        %%
        newBinEdges(1)   = newBinEdges(1)-1e-5;
        newBinEdges(end) = newBinEdges(end)+1e-5;        
        nBins = length(newBinEdges)-1;
        cnt_binned = zeros(1,nBins);
        for bin_i = 1:nBins            
            idx =  dPhis > newBinEdges(bin_i) & dPhis <= newBinEdges(bin_i+1) ;
            cnt_binned(bin_i) = sum(dPhiNull_Count(idx));
        end
        assert( abs( sum(cnt_binned) - sum(dPhiNull_Count)) < 1e-5 );
    end
