function [U, Umax] = gravitationalPotential(M, wrapDims, pow )
    if nargin < 2
        wrapDims = [1];  % first dimension (orientation) is cyclic.
    end
    if nargin < 3
        pow = 1;
    end
%     function s = wrapDist(i,j, L)
%         s = abs(i-j);
%         s = min(s, L-s);        
%     end
%     distance = @(x1, x2) norm(x1-x2);

    U = 0;
    Umax = 0;
    sizeM = size(M);
    dim1 = size(M,1);
    
    inds_ND = findV(M);
    inds_1D = find(M(:));
    N = length(inds_1D);
    for i = 1:N
        for j = 1:i-1
            ind_i_ND = inds_ND(i,:);
            ind_j_ND = inds_ND(j,:);
            
            ind_i_1D = inds_1D(i);
            ind_j_1D = inds_1D(j);
            Mi = M(ind_i_1D);
            Mj = M(ind_j_1D);
            distVec = ind_i_ND - ind_j_ND;
%             for d = wrapDims
                s = abs(ind_i_ND(1) - ind_j_ND(1));
                distVec(1) = min(s, dim1-s);
%             end
            U = U - Mi * Mj / (mag( distVec ) ^ pow);            
            
            Umax = Umax - Mi * Mj;
        end
        
        
    end

    
end