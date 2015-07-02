function testGetBestSpikeOffsets
    hnd = dbOpenExpDb;
    S = load('cellsGroups_grating');
    gratingGroups = S.gratingGroups;
    
    GrpId = 192;
    cellId = 2;
    groupData = gratingGroups(findInStructArray(gratingGroups, 'Gid', GrpId));

    delay_ms = 30;
    Gid = groupData.Gid;
    Did = groupData.Did;
    fieldsToGet = {'DBL_TEMP_PERIOD_FRM', 'LNG_START_TICK', 'LNG_END_TICK' };
    [period_frm, startTick, endTick] = getFieldsFromDatabaseTable(hnd, fieldsToGet, 'TBL_GRATING_PRES', {'DATAFILE_ID', Did});

    period_ms = dbConvertTimeMeasures(Did, period_frm, 'frame', 'ms');
    [start_ms end_ms] = elements( dbConvertTimeMeasures(Did, [startTick, endTick], 'tick', 'ms') );

    nPres = ceil( dbConvertTimeMeasures(Did, [endTick-startTick], 'tick', 'ms') / period_ms );  
               
%         cellSpikes = cell(1, length(cellIds));

    spikeTimes_ms = dbGetSpikes(Gid, cellId, 'ms') - delay_ms;
    spikeTimes_ms = elementsInRange(spikeTimes_ms, [start_ms, end_ms]);
    

    syncTimes_ms = dbGetSyncs('Did', Did, 'ms');
    spikePosInFrames = binarySearch(syncTimes_ms, spikeTimes_ms, -1, 0.5);
    spikePosInPres   = binarySearch([1:nPres]*period_frm, spikePosInFrames, -1, 0.5); % --> which Presentation the spikes are in 
    
    
    figure(150); clf;
    subplot(3,2,[1 3]);
    plotSpikeRaster(spikePosInPres, 1, period_ms, 'b');
    
    [offsets, offsetSpikeTimes] = getBestSpikeOffsets(spikePosInPres, 1, 1/100);
    
    subplot(3,2,[2 4]);
    plotSpikeRaster(spikePosInPres, 1, period_ms, 'r', offsets);
    
    
    subplot(3,2,[5 6]);     % histogram
%     d = getBestBinSize(spikePosInPres, 1, 'bar');
%     d = 1/30;
%     N_bins = round( 1/d );

    msPerSecond = 1000;   periodsPerSec = period_ms/msPerSecond;
    [x_bins_cent, N_hist] = calcPSTH( mod(offsetSpikeTimes, 1), 1, nPres, periodsPerSec);
    
    [isSimple, F1oDC] = isSimpleCell(x_bins_cent, N_hist, 1);
    cellType = iff(isSimple, 'Simple', 'Complex');

            
    plotThisPSTH(x_bins_cent, N_hist, 1);
    title(['Group ' num2str(Gid), '.  Cell ' num2str(cellId) ': ' cellType ' cell (F1/DC = ' num2str(F1oDC) ')']);
            
%     [offsets, offsetSpikeTimes] = getBestSpikeOffsets(spikePosInPres, 1, 1/100, 'flag');
    

end


%     singleInds = findInStructArray(gratingGroups, 'gratingType', 'Single Grating', @strcmp);
%     flashInds = findInStructArray(gratingGroups, 'gratingType', 'Flashed Grating Batch', @strcmp);
%     orientationInds = findInStructArray(gratingGroups, 'gratingType', 'Orientation Batch', @strcmp);
%     spatialFreqInds = findInStructArray(gratingGroups, 'gratingType', 'Spatial Frequency Batch', @strcmp);
%     temporalFreqInds = findInStructArray(gratingGroups, 'gratingType', 'Temporal Frequency Batch', @strcmp);
%     freeInds = findInStructArray(gratingGroups, 'gratingType', 'Free Grating Batch', @strcmp);
%     fieldsToGet = 'DBL_ORIENTATION_DEGR', 'DBL_TEMP_PHASE_DEGR', 'DBL_TEMP_PERIOD_FRM', 'DBL_SPATIAL_PHASE_DEGR', 'DBL_SPATIAL_PERIOD_PIX', 


%   unique(getFieldsFromDatabaseTable(hnd, 'DBL_TEMP_PERIOD_FRM', 'TBL_GRATING_PRES')) =
%     3.7500
%     5.0000
%     6.0000
%     6.6667
%     7.5000
%    10.0000
%    12.0000
%    15.0000
%    20.0000
%    30.0000
%    40.0000
%    43.7370
%    48.0000
%    60.0000
%   120.0000
%   240.0000
%   600.0000


%             [ts, fs] = getSpikeDensity(spikeTimesInRange_ms, mspf_approx/2);
%             [freqs, amps] = findStrongestFrequencies(ts,fs-mean(fs), [], 1);
% 
%             [tmp, ind_freq] = max(amps .* freqs);   % (ideally, use maxfreqs alone.)
%             period_ms = 1/ind_freq;

%             for ai = 1:100
%                 spk = .05*randn(round(10+rand*10), 1) + .5;
%                 spikePosInPres1{ai} = (1.0)*(ai-1) + spk;
%                 spikePosInPres2{ai} = (1.01)*(ai-1) + spk;
%             end
%             spikePosInPres1 = sort( [ vertcat(spikePosInPres1{:}) ] );
%             spikePosInPres2 = sort( [ vertcat(spikePosInPres2{:}) ] );
            

% success stories
% 6(103)  8(105) 9(106)






%         for ai = 1:100
%             spk = .1*randn(round(1+rand*5), 1) + .5;
%             spikePosInPres{ai} = (1.0)*(ai-1) + spk;
%             spikePosInPresO{ai} = (1.01)*(ai-1) + spk;
%         end
%         spikePosInPres  = sort( [ vertcat(spikePosInPres{:}) ] );
%         spikePosInPresO = sort( [ vertcat(spikePosInPresO{:}) ] );    
%     
%     
%     figure(154); clf;
%     subplot(2,2,1);
%     plotSpikeRaster(spikePosInPres, 1, period_ms, 'b');
% 
%     subplot(2,2,3);
%     plotSpikeRaster(spikePosInPresO, 1, period_ms, 'b');
%     
%     [offsets, offsetSpikeTimes] = getBestSpikeOffsets(spikePosInPresO, 1, 1/100);
%     
%     subplot(2,2,2);
%     plotSpikeRaster(spikePosInPresO, 1, period_ms, 'r', offsets);
%     
%     
%     subplot(2,2,4);     % histogram
% %     d = getBestBinSize(spikePosInPres, 1, 'bar');
%     d = 1/30;
%     N_bins = round( 1/d );

%     msPerSecond = 1000;   periodsPerSec = period_ms/msPerSecond;
%     [x_bins_cent, N_hist] = calcPSTH( mod(offsetSpikeTimes, 1), N_bins, periodsPerSec, nPres);
%     [isSimple, F1oDC] = isSimpleCell(x_bins_cent, N_hist, 1);
%     cellType = iff(isSimple, 'Simple', 'Complex');
% 
%             
%     plotThisPSTH(x_bins_cent, N_hist, 1);
%     title(['Group ' num2str(Gid), '.  Cell ' num2str(cellId) ': ' cellType ' cell (F1/DC = ' num2str(F1oDC) ')']);
%             
%     [offsets, offsetSpikeTimes] = getBestSpikeOffsets(spikePosInPresO, 1, 1/100);
