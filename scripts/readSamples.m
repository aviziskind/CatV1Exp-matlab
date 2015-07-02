function samples = readSamples(file_info, channelIds, sampleStart, nSamplesRead, skip_flag)
    % file_info should have 4 fields : 'fid', 'nChannelsTot', 'bytesPerValue', 'dataPrecision'.
    useReadSkip = exist('skip_flag', 'var') && ~isempty(skip_flag);

    fid = file_info.fid;
    nChannelsTot = file_info.nChannelsTot;  
    bytesPerValue = file_info.bytesPerValue;    
    
    nChannelsRead = length(channelIds);
    if useReadSkip                
        assert( (length(channelIds) == 1) ); % can only handle 1 channel   (( || all(diff(channelIds) == 1) ); % 1 channel or consecutive channels ))
        firstChannelId = channelIds(1);
        nBytesToSkip = (nChannelsTot-1)*bytesPerValue;
        startByte = (sampleStart-1)*nChannelsTot*bytesPerValue + (firstChannelId-1)*bytesPerValue;        
        fseek(fid, startByte, 'bof');  error(ferror(fid)); 
        samples = fread(fid, [nChannelsRead nSamplesRead], file_info.dataPrecision, nBytesToSkip);  

    else        
        startByte = (sampleStart-1)*nChannelsTot*bytesPerValue;
        fseek(fid, startByte, 'bof');  error(ferror(fid));
        samples = fread(fid, [nChannelsTot, nSamplesRead], file_info.dataPrecision);  error(ferror(fid));
        samples = samples(channelIds,:);
    end       

end
