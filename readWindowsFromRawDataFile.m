function [data, dfStruct] = readWindowsFromRawDataFile(df, t_starts, t_length, nChannels, channelIds)
    % data is arranged by [t_start..t_end] x [channel_1 .. channel_n] x [t_wind_1 .. t_wind_n] 

    keepOpen = false;
    if isstruct(df)
        dfname = df.dfname;
        nChannels = df.nChannels;
        channelIds = df.channelIds;
        fid = df.fid;
        keepOpen = true;
                
    elseif ischar(df)
        dfname = df;
            
        if nargin < 4  % lookup number of channels in database
            nChannels = getFieldsFromDatabaseTable([], 'LNG_N_CHANNELS', 'TBL_DATA_FILES', {'TXT_DATAFILE_NAME', dfname});
        end        
        if nargin < 5
            % find groupId, then find which channels.
            hnd = dbOpenExpDb;
            Did = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES', {'TXT_DATAFILE_NAME', dfname});
            Gid = dbLookup('Gid',  'Did', Did);
            if length(Gid) > 1
                error('There are multiple Groups associated with this datafile. Please specify which channels you want manually');
            end
            nChannelFields = 16;
            channelsUsedFieldnames = arrayfun(@(ch_id) ['BLN_CH_' num2str(ch_id, '%02d')], 1:nChannelFields, 'un', 0);
            [channelsUsed{1:nChannelFields}] = getFieldsFromDatabaseTable(hnd, channelsUsedFieldnames, 'TBL_GROUPS', {'GROUP_ID', Gid});
            channelIds = find([channelsUsed{:}]);
            warning('This function will be faster if you specify which channels you want');
        end
        if (length(df) < 3) || ~strcmp(df(2:3), ':\') % ie, doesn't contain full path.
            pathname = rawDataDir(df);
            df = [pathname df];
        end        
        fid = fopen(df);
    end

    if nargin < 2
        t_starts = 1;
        t_length = 1000;
    end
        
    nWindows = length(t_starts);

    dataPrecision = 'int16';
    bytesPerValue = 2;    
    
    data = zeros(nChannels, t_length, nWindows, dataPrecision); % this is the order in which the data is stored (first channel order, then time order).
    for tw_i = 1:nWindows
        fseek(fid, (t_starts(tw_i)-1)*nChannels*bytesPerValue, 'bof');
        error(ferror(fid))
        data(:,:,tw_i) = fread(fid, double([nChannels, t_length]), dataPrecision);
    end    
        
    data = data(channelIds,:,:);  % discard unwanted channel data
    data = permute(data, [2 1 3]); % put voltage trace from each channel in a continuous format.

    if nargout > 1 
        dfStruct = struct('dfname', dfname, 'nChannels', nChannels, 'channelIds', channelIds, 'fid', fid);
        keepOpen = true;
    end

    if ~keepOpen
        fclose(fid);
    end

end

%{
    % this takes longer
    tic
    data2 = zeros(4, t_length, nWindows);
    for tw_i = 1:nWindows
        fseek(fid, (t_starts(tw_i)-1)*nChannels*bytesPerValue, 'bof');
        error(ferror(fid))
        A = fread(fid, [4, t_length], '4*int16', bytesPerValue);
        data2(:,:,tw_i) = A;
    end



data = readRawDataFile(dfname, [1, 501, 1001], 1000);
figure(1); clf; plot(1:1000, data(:,:,1) ); xlim([0 2000])
figure(2); clf; plot(501:1500, data(:,:,2) ); xlim([0 2000])
figure(3); clf; plot(1001:2000, data(:,:,3) ); xlim([0 2000])

figure(1); clf; plot(1:1000, data(:,:,1), 'o' ); 
hold on
plot(501:1500, data(:,:,2), 's' ); 
figure(3); clf; plot(1001:2000, data(:,:,3) ); xlim([0 2000])

%}





