function mid_getCellNonLinearity(Gid, cellId)

    s = calculatePSTH_STAs_OSP_ForOneCell(Gid, cellId);
    OSP = s.OSP;

    mid_fileName = mid_getPreferredMIDfile(Gid, cellId);
    S_mid = load(mid_fileName);
    
    for i 







end