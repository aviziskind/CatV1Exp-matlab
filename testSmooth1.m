% testSmooth1

A = rand(36,10);

w = .8;
circFLAG = 1;
omitFLAG = 1;
% B12 = gaussSmooth(gaussSmooth(A, w, 1), w, 2, circFLAG);
% B21 = gaussSmooth(gaussSmooth(A, w, 2, circFLAG), w, 1);

B12 = gaussSmooth(gaussSmooth(A, w, 1), w, 2, circFLAG);
B21 = gaussSmooth(gaussSmooth(A, w, 1, [], omitFLAG), w, 2, circFLAG, omitFLAG);

figure(1);
clf;
subplot(1,3,1); imagesc(A);
subplot(1,3,2); imagesc(B12);
subplot(1,3,3); imagesc(B21);
d = max(B12(:)-B21(:));
title(sprintf('%.2g', d))


% Y2 = gaussSmooth(Y1, w, dim, circularFlag)