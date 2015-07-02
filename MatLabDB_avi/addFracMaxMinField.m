function pairList = addFracMaxMinField(pairList, allOSPs)
    % helper function for generateGratingDatafiles.

    % for each pair: compute the fraction of each cell's response at
    % the max ori/sp of the other cell vs its own max.        
    for i = 1:size(pairList,1);
        cell1 = pairList(i).cell1;
        cell2 = pairList(i).cell2;            
        mR1 = mean(allOSPs(cell1).R, 3);  %*average* across phases!
        mR2 = mean(allOSPs(cell2).R, 3);
        
        fR1 = mR1/max(mR1(:));
        fR2 = mR2/max(mR2(:));
        fMin = min(fR1, fR2);
        fMaxMin = maxElement(fMin);
        pairList(i).fracMaxMin = fMaxMin;
        
    end
        
end

%         [m1, indmax1] = maxElement( mR1 ); [ori1, sp1] = elements(indmax1);
%         [m2, indmax2] = maxElement( mR2 ); [ori2, sp2] = elements(indmax2);
%         pairList(i).frac1AtMax2 = mR1(ori2, sp2)/m1;
%         pairList(i).frac2AtMax1 = mR2(ori1, sp1)/m2;
