function r = decorrelatePSTHs_helper(R, nBinsPerFrame, r_framesConsidered, A)

    [nBinsPerExtFrame, nStims] = size(R);
    nFramesPerExtFrame = nBinsPerExtFrame/nBinsPerFrame;
    
    R_framesConsidered = 1:nFramesPerExtFrame;
    R_binsConsidered = crossOp([1:nBinsPerFrame]', @plus, [R_framesConsidered-1]*nBinsPerFrame );    
    r_binsConsidered = crossOp([1:nBinsPerFrame]', @plus, [r_framesConsidered-1]*nBinsPerFrame );
    
    R_shapedForMtx = reshape(R(R_binsConsidered(:),:), [nBinsPerFrame, length(R_framesConsidered)*nStims])';

    % Do the actual decorrlation:
%     disp('linsolve: (beginning linear solving ... wish me luck!)');
%     pause(1);
%     tic;
%     r_shapedForMtx2 = linsolve(A, R_shapedForMtx);
%     toc;

%     disp('mldivide: (beginning matrix inversion ... wish me luck!)');
%     pause(1);
%     tic;
    r_shapedForMtx = A \ R_shapedForMtx;
%     toc;    
    
    r  = reshape(r_shapedForMtx',  [nBinsPerFrame*length(r_framesConsidered), nStims]);
%     r2 = reshape(r_shapedForMtx2', [nBinsPerFrame*length(r_framesConsidered), nStims]);
%     save solvedMatrices_inv.mat r r2;
    
    % optionally: put r back into matrix the same size as the original R:
    r_big = zeros(size(R));
    r_big(r_binsConsidered(:),:) = r;
    r = r_big;

end



%     if ~exist('eigThreshold', 'var') || (eigThreshold == 0)
%         invA = inv(A);
%     else
%         [V,D] = eig(A);
%         invA = (V * invDiagThreshold(D, eigThreshold) * inv(V));
%     end


