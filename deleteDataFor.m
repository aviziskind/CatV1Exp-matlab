function deleteDataFor(Gids)

    whichSortingFeatures = 'all'; 
%     whichSortingFeatures = 'current'; 

    df_sharpWaveFileTypes = {'sharpWaves'};

    df_ECprocessingFileTypes = {...
        'meansCov',      ...         
        'startEnd',      ... 
        'properties',    ... 
        'waveforms',     ... 
        'waveforms_ext', ... 
        'waveformsPCA',  ... 
        'waveformsGLF'};

    df_ICprocessingFileTypes = {...
        'meanVarIC',     ...
        'startEndIC',    ... 
        'propertiesIC',  ... 
        'waveformsIC',   ... 
    };
    
    
    ECspikeProcessing_fileTypes = {...        
        'kNN',           ... 
        'clusters',      ... 
        'clusterData1',  ... 
        'clusterPrunings', ...
        'prunedClustersStats', ...
        'clustersPruned', ...
        'clusterData2',  ... 
        'clusterMerges', ...         
        'mwvfm_clust',   ... 
        'mwvfm_clustP',  ... 
        'mwvfm_cell', ...
        'cells', ...
        'IC_EC_match',  ...
    };

    spikePruning_fileTypes = {...        
        'clusterData1',  ... 
        'clusterPrunings', ...
        'prunedClustersStats', ...
        'clustersPruned'};
    
    
    sortingSpecific_fileTypes = setdiff(ECspikeProcessing_fileTypes, 'kNN');
    
%         allFileTypes = [df_sharpWaveFileTypes, df_ECprocessingFileTypes, df_ICprocessingFileTypes, ECspikeProcessing_fileTypes];

%     allFileTypes = [ECspikeProcessing_fileTypes];
%     allFileTypes = [df_ECprocessingFileTypes, ECspikeProcessing_fileTypes];
%     allFileTypes = [df_sharpWaveFileTypes, df_ICprocessingFileTypes];
%     allFileTypes = [df_ECprocessingFileTypes, ECspikeProcessing_fileTypes];
        
    allFileTypes = [spikePruning_fileTypes];
    

    
    switch whichSortingFeatures
        case 'all', 
            features_list_S = dir('C:\ExperimentDB\Spikes\clusters');
            features_list = {features_list_S.name};
            fet_strs = features_list( cellfun(@(s) isstrprop(s(1), 'alpha'), features_list ) );
        case 'current',
            fet_strs = curSortingFeatures('');
    end
    
    if ischar(fet_strs)
        fet_strs = {fet_strs};
    end
    
        
    progressBar('init-', length(Gids));
    for i = 1:length(Gids)
        Gid = Gids(i);
        
        for j = 1:length(allFileTypes)
            fileType = allFileTypes{j};
            if any(strcmp(sortingSpecific_fileTypes, fileType))
                
                for fet_i = 1:length(fet_strs)                                    
                    fn = getFileName(fileType, Gid, [], struct('sortingFeatures', fet_strs{fet_i}) );

                    if exist(fn, 'file');                
                        delete(fn);
                    end                                
                end        
            else
                fn = getFileName(fileType, Gid);
            
                if exist(fn, 'file');                
                    delete(fn);
                end            
            end
        end
        progressBar;
    end
        





end

%{

                if k == 0
                    fprintf('Matlab : '); tic;
                    delete(fn); toc;
%                     k = 1;
                elseif k == 1
                    fprintf('Dos : '); tic;
                    system(['del ' fn]); toc;
                    k = 0;
                end                    

%}