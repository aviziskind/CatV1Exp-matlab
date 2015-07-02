function viewExperimentChronology

%     plotMode = 'grating'; %options: 'all', 'grating', 'noise';
    plotMode = 'all'; %options: 'all', 'grating', 'noise';
    dateMask = 'mm/dd/yyyy HH:MM:SS PM';

    hnd = dbOpenExpDb;
    [Dids, stimulusIds, dates_did] = getFieldsFromDatabaseTable(hnd, {'DATAFILE_ID', 'STIMULUS_TYPE_ID', 'DTM_CREATED'}, 'TBL_DATA_FILES');
    N = length(Dids);
    getStimName = @(id) getFieldsFromDatabaseTable(hnd, 'TXT_STIMULUS_TYPE', 'TBL_STIMULUS_TYPES', {'STIMULUS_TYPE_ID', id});
    stimTypes = {'grating', 'noise',  'mseq', 'movie' };
    stimTypeIds = {[1,2,3,6,7,9], 4, 5, [8, 10]};

    gratingStimIds = [1,2,3,6,7,9];
    nGratingStims = arrayfun(@(id) nnz(stimulusIds == id), gratingStimIds);
    gratingStimLabels = arrayfun(@(id, n) [cell2mat(getStimName(id)) sprintf(' (%d)', n )], gratingStimIds, nGratingStims, 'un', 0);

    noiseStimLabels = {'Noise'};
    mseqStimLabels = {'Mseq'};
    movieStimLabels = {'Movie', 'Movie2'};
    allStimLabels = {gratingStimLabels, noiseStimLabels, mseqStimLabels, movieStimLabels};
    
    [stimIdTables, presNames] = deal( cell(1,10) );
    for si = 1:length(stimTypes)
        stimIdTables(stimTypeIds{si}) = {['TBL_' upper(stimTypes{si}) '_PRES']};
        presNames(stimTypeIds{si}) = {[upper(stimTypes{si}) '_PRES_ID']};
    end
    
    if ~exist('exp_chron.mat', 'file')

        dates_dids = datenum(dates_did, dateMask);    

    %     [presIds, dates_pres] = deal( zeros(1,N) );
        [dates_pres1, dates_presN] = deal( zeros(1,N) );
        progressBar('init-', N, 50);
        for di = 1:N
            progressBar(di);
            si = stimulusIds(di);
            [presIds, dates_pres] = getFieldsFromDatabaseTable(hnd, {presNames{si}, 'DTM_CREATED'}, stimIdTables{si}, {'DATAFILE_ID', Dids(di)}, 'LNG_PRESENT_NO');
            [tmp, id_min] = min(presIds);
            [tmp, id_max] = max(presIds);
            dates_pres1(di) = datenum(dates_pres(id_min), dateMask);
            dates_presN(di) = datenum(dates_pres(id_max), dateMask);
        end
        progressBar('done');
        
        save('exp_chron.mat', 'dates_*');
    else
        load('exp_chron.mat')
    end
    3;
%     dates = getFieldsFromDatabaseTable(hnd, {'STIMULUS_TYPE_ID', 'DTM_CREATED'}, 'TBL_DATA_FILES', {'DATA_FILE_ID', 1252});
    
    if strcmp(plotMode, 'all')
        figure(15); clf; hold on;
        nStimTot = 10;
        stim_i = 1;        
        ytkLabels = cell(1,nStimTot);
        for s_type_i = 1:length(stimTypeIds)
            stimTypeIds_i = stimTypeIds{s_type_i};
            for s_subtype_i = 1:length(stimTypeIds_i)                        
        
                idx = stimulusIds == stimTypeIds_i(s_subtype_i);
                plot(dates_pres1(idx), ones(nnz(idx))*stim_i, [color_s(stim_i) 'o']);
                plot(dates_presN(idx), ones(nnz(idx))*stim_i, [color_s(stim_i) 's']);
        %         plot(dates_dids(idx), ones(nnz(idx))*stim_i, [color_s(stim_i) '*']);
                line([dates_pres1(idx); dates_presN(idx)], stim_i*[1;1], 'color', color_s(stim_i));        
                ytkLabels{stim_i} = allStimLabels{s_type_i}{s_subtype_i};
                stim_i = stim_i+1;
                
            end
        end
        set(gca, 'ytick', 1:nStimTot, 'ytickLabel', ytkLabels);
        ylim([0, nStimTot+1])
        datetick
        
    elseif strcmp(plotMode, 'grating')
        figure(15); clf; hold on;
        for gi = 1:length(gratingStimIds)
            idx = stimulusIds == gratingStimIds(gi);
            plot(dates_pres1(idx), ones(nnz(idx))*gi, [color_s(gi) 'o']);
            plot(dates_presN(idx), ones(nnz(idx))*gi, [color_s(gi) 's']);
    %         plot(dates_dids(idx), ones(nnz(idx))*gi, [color_s(i) '*']);
            line([dates_pres1(idx); dates_presN(idx)], gi*[1;1], 'color', color_s(gi));        
        end
        set(gca, 'ytick', 1:length(gratingStimIds), 'ytickLabel', gratingStimLabels);
        ylim([0, length(gratingStimIds)+1])
        datetick
        
    end

    n20min = datenum(0,0,0,0,20,0);
    
end
%     for i = 1:
    
    
    
    
        
%     low = min(datenums);
%     high = max(datenums);
%     lowVec  = datevec(low);    lowYear = lowVec(1);
%     highVec = datevec(high);   highYear = highVec(1)+1;
%     
%     yearTicks = lowYear:highYear;
%     tickVals = zeros(size(yearTicks));
%     for i = 1:length(yearTicks)
%         tickVals(i) = datenum([yearTicks(i) 0 0 0 0 0]);
%     end
    

%     stimulusIds2 = zeros(size(stimulusIds));
% 
%     stimTypes = {'grating', 'noise',  'mseq', 'movie' };
%     stimTypeIds = {[1,2,3,6,7,9], 4, 5, [8, 10]};
%     for i = 1:4
%         for j = 1:length(stimTypeIds{i})
%             stimulusIds2(stimulusIds == stimTypeIds{i}(j)) = i;
%         end
%     end
        
%     colOrd = zeros(1,10);
    
%     colorOrder = color(colOrd);
%     h = plot(datenums, stimulusIds2, 'o');
%         
%     datetick;
% %     set(gca, 'xtick', tickVals, 'xticklabel', yearTicks);
%     set(gca, 'ytick', 1:4, 'yticklabel', stimTypes); 
    
% end

% 2669