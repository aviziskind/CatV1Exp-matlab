function dbDoAllDBCorrections

    % small corrections
    dbCorrectFrameLengths; % two lng_frames_per_update  fields corrected in noise presentation table    
    dbUpdateMoviePathNames;    
    dbCorrectGratingStaticFields; % some drifting gratings have static==true, some flashed have static==false.
    dbFixMovieTableFields; % some N_Frames fields are incorrect.
    dbFixTempFreqBatchStimIds % some StimulusTypeIds are mislabeled (i think?).
    dbFixElectrodeTypeId;
    
    % major corrections
    dbInsertMissingRecords; % insert 2 missing records into movie_pres table           
    dbFixSustainedDisplayedFrms;
    dbRedoAllTbTe;    
    
    % pruning bad records
    dbRemoveGroupsWithNoSpikes;
    dbRemoveDiodeTests;
    dbRemoveOrphanedGidsDids;  
%     dbRemoveEntriesWithNoSyncs;    
    


end