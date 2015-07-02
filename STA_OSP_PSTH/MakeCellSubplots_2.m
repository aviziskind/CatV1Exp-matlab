% benchmark for my code: should have same output as MakeCellSubplots, but uses my functions
% instead of Sergei's
tic;

Gid = 4470;
cellId = 2;

% % Gid = 4712;
% % Gid = 4726;

comp_ORIGINAL = 0;
comp_IMPROVED = 1;
computation = comp_ORIGINAL;
% computation = comp_IMPROVED;

BkgrSkip_ms = 200;
timeWindow_ms = [0];

[spkTsRelToFrame_ms, bckgRate] = getParsedSpikes('timing', Gid, cellId, BkgrSkip_ms);    

Mids = dbLookup('Mid',  'Gid', Gid);


%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% PSTH %%%%%%%%%%%%%%

nTotalFrames = length( spkTsRelToFrame_ms );
frameLength_ms = getFrameLength('Gid', Gid);


[PSTH_bins, PSTH_vals] = calcPSTH( spkTsRelToFrame_ms, frameLength_ms, nTotalFrames);

if (computation == comp_ORIGINAL) 
    nSpikesEachFrame = getParsedSpikes('frame', Gid, cellId);
elseif (computation == comp_IMPROVED)
    [timeWindow, windowProfile] = getBestTimeWindowFromPSTH(PSTH_bins, PSTH_vals, 'newPlot', bckgRate);
    relContrOfFrameToSpike = getParsedSpikes('frame', Gid, cellId, timeWindow, windowProfile );
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% OSP %%%%%%%%%%%%%%

if (computation == comp_ORIGINAL) 
    [OSP, R_full, oris, sps, phs] = getOriSpfPhaseProfile_simple(Gid, [nSpikesEachFrame{:}]);
elseif (computation == comp_IMPROVED)
    [OSP, R_full, oris, sps, phs] = getOriSpfPhaseProfile_simple(Gid, relContrOfFrameToSpike);
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% STA %%%%%%%%%%%%%%


if (computation == comp_ORIGINAL) 
    timeWindow = {[0 1], 'frame'};
    STAs = getSTAforCell(Gid, cellId, timeWindow, [], [nSpikesEachFrame{:}] );
elseif (computation == comp_IMPROVED)
    STAs = getSTAforCell(Gid, cellId, timeWindow, windowProfile, [relContrOfFrameToSpike{:}]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOT EVERYTHING %%%%%%%

figure(2);
ha_0 = subplot(1,3,1);
plotThisPSTH(PSTH_bins, PSTH_vals, frameLength_ms, bckgRate);
axis square;

subplot(1,3,2);
imagesc(oris, sps, sum(OSP, 3));
axis square;

subplot(1,3,3);
imagesc(STAs(:,:,1)); 
axis square xy;


hnd = dbOpenExpDb;

fieldnames = {'TBL_ELECTRODES.ELECTRODE_TYPE_ID', 'TBL_ANIMALS.TXT_LAB_NAME', 'TBL_PENETRATIONS.PENETRATION_ID', 'TBL_LOCATIONS.LOCATION_ID', 'TBL_DATA_FILES.LNG_DATAFILE_NO',  'TBL_DATA_FILES.DBL_SAMPLING_RATE_HZ'};

T1 = {'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'};               
T2 = {T1, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};               
T3 = {'TBL_ANIMALS', 'ANIMAL_ID', T2, 'TBL_ELECTRODES'};
T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'};
T5 = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
T6 = {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  T5, 'TBL_LOCS_FILES_LINKS'};
T7 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', T6, 'TBL_LOC_DEPTHS'};
T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};
T9 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', T8, 'TBL_ELECTRODES'};
joinedTables = T9;

% joinedTables = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', {{{{'TBL_ANIMALS', 'ANIMAL_ID', {{'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'}, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'}, 'TBL_ELECTRODES'}, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'}, 'TBL_LOCATIONS', 'LOCATION_ID', {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'}, 'TBL_LOCS_FILES_LINKS'}, 'TBL_LOC_DEPTHS'}, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'}, 'TBL_ELECTRODES'};


% criterea = {'TBL_GROUPS.GROUP_ID', Gid};
criterea = {'TBL_ELECTRODES.ELECTRODE_TYPE_ID', 2; 'TBL_GROUPS.GROUP_ID', Gid};
[elecType, CatName, PenetrID, LocID, FileNmb, SampRateHz] = getFieldsFromDatabaseTable(hnd, fieldnames, joinedTables, criterea);

str_1 = [ CatName '   PenID: ' num2str(PenetrID) '   LocID: ' num2str(LocID) '   Gid: ' num2str(Gid)  ];
str_2 = ['File #: ' num2str(FileNmb) '  Cell# ' num2str(cellId) '  N spikes: ' num2str(length(spkTsRelToFrame_ms)) ];

suptitle_2([str_1 sprintf('\n') str_2]);

toc;       
     