% testSTCMdecorrelation

Gid = 2078;
cellId = 1;
% $STA = getSTAforCell(Gid, cellId, {30, 'ms'}, {30, 'ms'});RS
rf = STCM_decorrelation(Gid, cellId);
figure; imagesc(rf);

