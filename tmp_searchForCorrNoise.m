function tmp_searchForCorrNoise
    Gid = 4470;

    [M,C] = getChannelMeansAndCovariance(4470);
    [V,D] = eig(C); 
    EvectMax = V(:,end);
    invCovMtx = inv(C);

    siteData = siteDataFor('Gid', Gid);
    dfInfo = siteData.dataFileInfo;
    nChannelsTot = dfInfo.nChannels;
    nSamplesTot = dfInfo.filesize/(2*nChannelsTot);            

    channelIds = dfInfo.channelIds;
    nChannelsToUse = length(channelIds);
    
    chunkSize = 500;
    nChunks = floor(nSamplesTot/chunkSize)*2-1;
    
    dfname = [dfInfo.dataFileName 'F'];
    pathname = rawDataDir(dfname);    
    dfname = [pathname dfname];
    
    chunk_ccs = zeros(4,4,nChunks);
    chunk_max_cc = zeros(1,nChunks);
    id_max = zeros(2,nChunks);
    
    progressBar('init', nChunks, 30);
    for ci = 1:nChunks
        progressBar;
        t_starts = 1 + chunkSize*(ci-1)/2;  % blank
        
        chunk_i = readWindowsFromRawDataFile(dfname, t_starts, chunkSize, nChannelsTot, channelIds);
        chunk_i = double(chunk_i);
%         chunk_i = bsxfun(@minus, chunk_i, mean(chunk_i, 1));
        [aboveTh_tmp, md] = aboveEllipsoidalThreshold(chunk_i', 4, EvectMax, invCovMtx);
        if any(md > 4)
            continue;
        end
        
        cc = corr(chunk_i);
        chunk_ccs(:,:, ci) = cc;
        cc_l = tril(cc, -1);
        [chunk_max_cc(ci), id_max(:,ci)] = maxElement( cc_l );        
    end


    3;
    
end
