function dbCorrectGratingStaticFields
    hnd = dbOpenExpDb;

    
%   {'Single Grating', 'Flashed Grating Batch', 'Orientation Batch', 'Spatial Frequency Batch', 'Temporal Frequency Batch', 'Free Grating Batch'}
%    1, 3, 2, 6, 7, 9, 10

    driftingStimIds = [1, 2, 6, 7, 9, 10];
    staticStimIds = 3;
    
    driftingDids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES', {'STIMULUS_TYPE_ID', driftingStimIds});
    driftingDids = unique(driftingDids);
    staticDids = getFieldsFromDatabaseTable(hnd, 'DATAFILE_ID', 'TBL_DATA_FILES', {'STIMULUS_TYPE_ID', staticStimIds});    
    staticDids = unique(staticDids);
    
    % find any drifting grating presentations with STATIC = TRUE
    [wrongDriftingPresIds, wrongDriftingDids] = getFieldsFromDatabaseTable(hnd, {'GRATING_PRES_ID', 'DATAFILE_ID'}, 'TBL_GRATING_PRES', ...
        {'DATAFILE_ID', driftingDids; 'BLN_STATIC_GRATING', -1});    
    wrongDriftingDids = unique(wrongDriftingDids);
    
    fprintf('There were %d / %d DRIFTING grating experiments with these DatafileIds incorrectly have static == true : \n ', length(wrongDriftingDids), length(driftingDids));
    disp(wrongDriftingDids')    
    
    % find any flashed grating presentations with STATIC = FALSE
    [wrongStaticPresIds, wrongStaticDids] = getFieldsFromDatabaseTable(hnd, {'GRATING_PRES_ID', 'DATAFILE_ID'}, 'TBL_GRATING_PRES', ...
        {'DATAFILE_ID', staticDids; 'BLN_STATIC_GRATING', 0});
    wrongStaticDids = unique(wrongStaticDids);

    fprintf('\nThere were %d / %d FLASHED grating experiments with the these DatafileIds incorrectly have static == false : \n ', length(wrongStaticDids), length(staticDids));
    disp(wrongStaticDids')        % 71 / 262

    display('Correcting bad drifting grating fields');
    % correct wrong drifting gratings;
    progressBar('init-', length(wrongDriftingPresIds));
    for i = 1:length(wrongDriftingPresIds)
        progressBar(i);        
        updateValueInDatabaseTable(hnd, 0,  'BLN_STATIC_GRATING', 'TBL_GRATING_PRES', {'GRATING_PRES_ID',  wrongDriftingPresIds(i)});
    end
    progressBar('done');

    display('Correcting bad flashed grating fields');
    % correct wrong static gratings;    
    progressBar('init-', length(wrongStaticPresIds));
    for i = 1:length(wrongStaticPresIds)
        progressBar(i);
        updateValueInDatabaseTable(hnd, -1,   'BLN_STATIC_GRATING', 'TBL_GRATING_PRES', {'GRATING_PRES_ID',  wrongStaticPresIds(i)});
    end
    progressBar('done');
    
    checkTmpFields  = true;
    if checkTmpFields  % don't need these, because all frame periods are correct.    
        wrongDriftingPresIds2 = getFieldsFromDatabaseTable(hnd, 'GRATING_PRES_ID', 'TBL_GRATING_PRES', ...
            {'DATAFILE_ID', driftingDids; 'DBL_TEMP_PERIOD_FRM', 50000}); % empty
        wrongStaticPresIds2 = getFieldsFromDatabaseTable(hnd, 'GRATING_PRES_ID', 'TBL_GRATING_PRES', ...
            {'DATAFILE_ID', staticDids; 'DBL_TEMP_PERIOD_FRM', {'<', 50000} });  % empty

        if ~isempty(wrongDriftingPresIds2)    
            setFieldsInDatabaseTable(hnd, unknownTmp, 'DBL_TEMP_PERIOD_FRM', 'TBL_GRATING_PRES', {'GRATING_PRES_ID', wrongDriftingPresIds2});
        end
        if ~isempty(wrongStaticPresIds2)        
            setFieldsInDatabaseTable(hnd, 50000, 'DBL_TEMP_PERIOD_FRM', 'TBL_GRATING_PRES', {'GRATING_PRES_ID', wrongStaticPresIds2});
        end
    end

    
    
%     updateValueInDatabaseTable(hnd, '9/25/2002 00:00:01', 'DTM_CREATED', 'TBL_GRATING_PRES', {'GRATING_PRES_ID', 166228}, 'DATE')
end



% wrongDriftingDids
%   [374  383  392 1309 1466 1482 1485 1486 1487 1777 3230 3231 3232 3233]
% wrongStaticDids
%   [425  437  438  444  445  452  453  455  461  462  466  467  471  472  475  476  477
%   478  487  488  489  492  493  494  498  503  504  505  506  510  511  512  513  517
%   519  523  524  526  527  535  536  543  544  545  546  548  549  556  557  558  559
%   560  561  565  566  567  572  575  576  579  580  581  585  586  587  588  616  617]


% wrongStaticGids = [ ...
%  481 463 432 465 466 473 474 476 483 484 488 490 492 493 494 495 496 ...
%  497 502 503 504 505 506 507 511 514 515 516 517 520 521 522 523 525 ...
%  526 528 530 529 531 537 538 543 544 545 546 453 454 549 550 551 552 ...
%  553 554 558 559 560 563 565 566 567 568 569 573 574 575 576 598 599 ];
% wrongDriftingGids = [ ...
%  [337 264 282 1048 1141 338 265 283 1254 1265 1266 1267 1502 4177 4178 4179 ...
%   4180 4181 4182 4183 4184];
  
  