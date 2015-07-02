function [dphis, nullP] = deltaPhiNull(nPh, nSamples, allNPh, nSummed)
    % [dphis, nullP] = deltaPhiNull(nPhOptions, nSamplesOfEachNPh)

    if (nargin < 2) || isempty(nSamples)
        nSamples = ones(size(nPh));
    end
        
    if (nargin >= 3) && ~isempty(allNPh)
        additionalNPh = setdiff(allNPh, nPh);
        nPh = [nPh, additionalNPh];
        nSamples = [nSamples, zeros(size(additionalNPh))];
    end    
    
    if (nargin < 4)
        nSummed = ones(size(nPh));
    end

    
    if length(nPh) == 1

        d = 360/nPh;
        dphis = 0:d:180;
        nDph = length(dphis);

        nullP = ones(size(dphis));
        nullP(dphis == 0) = .5;
        nullP(dphis == 180) = .5;
        nullP = nullP/sum(nullP); % normalize 
        
        if nSummed == 2     
            %%
            pxp = nullP(:)*nullP(:)';            
            
            nDph_sing = nPh/2+1;
            nDph_double = nDph_sing*2-1;
                        
            idxs = -nDph_sing+1:nDph_sing-1;
            
            nullP2 = zeros(1, nDph_double);
            for i = 1:length(idxs)                
                nullP2(i) = sum( diag(pxp, idxs(i)) );                
            end            
            assert( sum(nullP2) == 1 );
            nullP = nullP2;
            dphis = 0:d/2:180;
        end
        
        nullP = nullP * nSamples;
        
    else
%         nSummed_arg = nSummed;
%         if length(nSummed_arg) < length(nPh)
%             nSummed_arg = nSummed*ones(size(nPh));
%         end
        
        if (nargin < 3) || isempty(allNPh)
            allNPh = nPh;
        end           
        
%         [all_dphis, all_P] = arrayfun(@deltaPhiNull, nPh, ones(size(nPh)), allNPh, nSummed_arg, 'un', 0);
        [all_dphis, all_P] = arrayfun(@(nph, nsmp) deltaPhiNull(nph, nsmp, nph, nSummed), nPh, ones(size(nPh)), 'un', 0);
        %%
%         for i = 1:length(nPh)
%             [all_dphis{i}, all_P{i}] = deltaPhiNull(nPh(i), 1, nPh(i), nSummed);
%         end
        %%
        
        dphis = unique([all_dphis{:}]);
        nullP = zeros(size(dphis));
        for i = 1:length(all_dphis)
            idx = arrayfun(@(ph) find(dphis == ph), all_dphis{i});
            nullP(idx) = nullP(idx) + all_P{i}*nSamples(i);
        end        
        
    end
        

end

%{
% 4 phases: single value
0    0.25
90   0.5
180  0.25

% 4 phases: two values
0    (0,0)            1 x (0.25)(0.25) = 1/16
45   (0,90; 90,0)     2 x (0.25)(0.5)  = 1/4
90   (0,180; 180,0;90,90)  3 ; 2*(0.25)(0.25)+(0.5)(0.5)  = 3/8
135  (90, 180;180,90) 1 x (0.25)(0.5)  = 1/4
180  (180, 180)       1 x (0.25)(0.25) = 1/16



% 8 phases: single value
0    0.125  1/8
45   0.25   1/4
90   0.25   1/4
135  0.25   1/4
180  0.125  1/8

% 8 phases: two values
0    (0,0)           1x (1/8)(1/8)  = 1/64
-22  (0,45;x2)       2x (1/8)(1/4)  = 1/16
45   (0,90;x2,45,45) 3 : 2*(1/8)(1/4) + (1/4)(1/4)  = 2/32 + 1/16 = 1/8
-67  (45,90;x2; 0,135;x2) 2x(1/4)(1/4) + 2*(1/8)(1/4) = 3/16
90   (0,180;x2, 45,135;x2, 90,90) = 2*(1/8)(1/8) + 2*(1/4)(1/4) + (1/4)(1/4) = 7/32
-112
135  
-157
180  0.5


%}