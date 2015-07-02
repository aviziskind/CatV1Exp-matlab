function checkSyncs(idType, idVal)
    Did = dbLookup('Did', idType, idVal);

    syncs = [0; dbGetSyncs('Did', Did, 'tick')];
    frameLength = getFrameLength('Did', Did, 'tick');
    
    dSyncs_frm = round ( diff(syncs) / frameLength );
    
    isFrame = double( (abs(dSyncs_frm-1) < .2) );
    isFrame(~isFrame) = dSyncs_frm(~isFrame);

    B = runLengths(isFrame);
    disp(B);



end