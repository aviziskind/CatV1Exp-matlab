function file_info = getOpenFileHandle(arg1, opt)
%   file_info = getOpenFileHandle(Gid, opt)
%   file_info = getOpenFileHandle(siteData, opt)
%   file_info = getOpenFileHandle(dataFileInfo, opt)
    filteredDatafilesUseOnlyRelevantChannels = 1;

    fileExt_default = '';
    outputPrecision_default = 'double';
    if (nargin == 2) && isfield(opt, 'fileExt')
        fileExt = opt.fileExt;
    else
        fileExt = fileExt_default;
    end
    
    if (nargin == 2) && isfield(opt, 'precision')
        outputPrecision = opt.precision;
    else
        outputPrecision = outputPrecision_default;
    end    
    
    if isnumeric(arg1)
        Gid = arg1;
        siteData = siteDataFor(Gid);
        dfInfo = siteData.dataFileInfo;
    elseif isstruct(arg1) && isfield(arg1, 'dataFileInfo')
        sd = arg1;
        dfInfo = sd.dataFileInfo;
    elseif isstruct(arg1) && isfield(arg1, 'dataFileName')
        dfInfo = arg1;
    end
        
    

    dfname = dfInfo.dataFileName;
	useLocalDir = ~isempty(fileExt);
    pathname = rawDataDir(dfname, useLocalDir);	
    
    dfname = [pathname dfname ];
    if ~isempty(fileExt)
        dfname = [dfname fileExt];
    end    
                       
    if ~isempty(strfind(dfname, 'datM'))
        src_precision = 'float32';
        bytesPerValue = 4;
        nChannelsTot_Orig = 1;
    else        
        src_precision = 'int16';
        bytesPerValue = 2;
        nChannelsTot_Orig = dfInfo.nChannels;
        if filteredDatafilesUseOnlyRelevantChannels  && ~isempty(fileExt)
            nChannelsTot = 4;
        else
            nChannelsTot = nChannelsTot_Orig;
        end        
    end
        
    dataPrecision = [ src_precision '=>' outputPrecision];
%     switch outputPrecision
%         case 'int16',  dataPrecision = 'int16';
%         case 'single', dataPrecision = 'int16=>single';
%         case 'double', dataPrecision = 'int16=>double';
%     end        
    
    if isfield(dfInfo, 'nSamples')
        nSamplesTot = dfInfo.nSamples;
    else
        nSamplesTot = dfInfo.filesize/(2*nChannelsTot_Orig*bytesPerValue);            
    end


    fid = fopen(dfname);  
    
    file_info.fid = fid;
    file_info.dfname = dfname;
    file_info.nChannelsTot = nChannelsTot;
    file_info.dataPrecision = dataPrecision;
    file_info.bytesPerValue = bytesPerValue;
    file_info.nSamples = nSamplesTot;
        
end
%     samples_IC_F = readSamples(file_info, channelId, sampleStart, bufferSize);