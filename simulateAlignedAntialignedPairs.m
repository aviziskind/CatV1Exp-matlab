function [nAA_exp, nAB_exp, nBB_exp] = simulateAlignedAntialignedPairs(nUnitsEachSite, probs)
    if nargin < 1
        probs = [0.71, 0.24];
    end
    %%
    maxNPerSite = max(nUnitsEachSite);
    possNUnits = 1:maxNPerSite;
%     C = arrayfun(@(n) nchoosek(1:n, 2), possNUnits, 'un', 0);
    %%
%     nPairsVNUnits = cellfun(@(x) size(x,1), C);
    [nAA_av_vs_nCellsInSite, nAB_av_vs_nCellsInSite, nBB_av_vs_nCellsInSite] = deal(zeros(1, length(possNUnits)));
    for i = 1:length(possNUnits)
        nUnits_i = possNUnits(i);
%         if nUnits_i == 1
%             continue;
%         end
                        
        possNaligned = 0:nUnits_i;
        possNanti_aligned = nUnits_i - possNaligned;
        
        probNaligned = binopdf(0:nUnits_i, nUnits_i, probs(1));
        [nAA, nAB, nBB] = getNumPairs(possNaligned, possNanti_aligned);
        %%
        nAA_av_vs_nCellsInSite(i) = sum(nAA .* probNaligned);
        nAB_av_vs_nCellsInSite(i) = sum(nAB .* probNaligned);
        nBB_av_vs_nCellsInSite(i) = sum(nBB .* probNaligned);
        
    end
    3;
    %%
    nAA_exp = sum( nAA_av_vs_nCellsInSite(nUnitsEachSite) );
    nAB_exp = sum( nAB_av_vs_nCellsInSite(nUnitsEachSite) );
    nBB_exp = sum( nBB_av_vs_nCellsInSite(nUnitsEachSite) );
    

end


function [nAA, nAB, nBB] = getNumPairs(possNaligned, possNanti_aligned)
    nAA = possNaligned .* (possNaligned-1) / 2;
    nAB = possNaligned .*possNanti_aligned;
    nBB = possNanti_aligned .* (possNanti_aligned-1) / 2;
end