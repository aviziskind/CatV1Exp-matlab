function [channelIds, nChannelsTot] = getChannelIds(arg1, voltageFilterMode)
    if ~strcmp(voltageFilterMode.filterName, 'none')
        channelIds = 1:4;
        nChannelsTot = 4;
    else
        if isnumeric(arg1)
            Gid = arg1;
            sd = siteDataFor(Gid);
            dfInfo = sd.dfInfo;
        elseif isstruct(arg1) && isfield(arg1, 'dataFileInfo');
            sd = arg1;
            dfInfo = sd.dataFileInfo;
        elseif isstruct(arg1) && isfield(arg1, 'channelIds');
            dfInfo = arg1;
        end            
        channelIds = dfInfo.channelIds;
        nChannelsTot = dfInfo.nChannelsTot;
    end
    
end