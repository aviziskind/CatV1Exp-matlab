function processDatafile_helper

    global CPU_id

     if exist([spikePath 'stop_running'], 'file')
        movefile([spikePath 'stop_running'], [spikePath 'stop_runnin'])
     end
   
     
%      allNs = [1:12, 14:2:20, 24:4:36];
%      allNs = [4:12, 14:2:20, 24:4:36];
%      allPCAcur = arrayfun(@(n) sprintf('PCAcur%d', n), allNs, 'un', 0);
%      allPCAcuw = arrayfun(@(n) sprintf('PCAcuw%d', n), allNs, 'un', 0);
%      allGLFcur = arrayfun(@(n) sprintf('GLFcur%d', n), allNs, 'un', 0);
%      allGLFcuw = arrayfun(@(n) sprintf('GLFcuw%d', n), allNs, 'un', 0);
%      
%      allFets = [allPCAcuw allPCAcur, allGLFcuw]; 
%      allFet_str = allFets;

%    allFets = [allGLFcur, allGLFcuw allPCAcur, allPCAcuw];
     
%      [allFets, allFet_str] = getAllSortingFeatures(4);

%         allFets = { ...
%         'Neg', 'GLFcuw4', 'GLFcuw8', 'GLFcuw12',    'GLFcuw16', ...
%         'GLFsnr3', ...
%          'PCAcur4', 'PCAcur8', 'PCAcur12', ...
%          'PCAcuw4', 'PCAcuw8', 'PCAcuw12' ...         
%          'PCAsnr3',  'PCAsur3', 'PCAsuw3', ...
%         };

        features_smoothN = { ...
            'Neg', 'PCAcuw2', 'PCAcuw3', 'PCAcuw4', 'PCAcuw5', 'PCAcuw6', 'PCAcuw8', 'PCAcuw10', 'PCAcuw12', 'PCAcuw16', ...
                   'PCAcur2', 'PCAcur3', 'PCAcur4', 'PCAcur5', 'PCAcur6', 'PCAcur8', 'PCAcur10', 'PCAcur12', 'PCAcur16', ...
                   'GLFcuw2', 'GLFcuw3', 'GLFcuw4', 'GLFcuw5', 'GLFcuw6', 'GLFcuw8', 'GLFcuw10', 'GLFcuw12', 'GLFcuw16', ...
                   'GLFcur2', 'GLFcur3', 'GLFcur4', 'GLFcur5', 'GLFcur6', 'GLFcur8', 'GLFcur10', 'GLFcur12', 'GLFcur16'};

        features_representative = { ...
         'Neg', 'PCAcuw4', 'PCAcuw8', 'PCAcuw12', 'PCAcuw16', ...
                'PCAcur4', 'PCAcur8', 'PCAcur12', 'PCAcur16', ...
                'PCAsur1', 'PCAsur2', 'PCAsur3', 'PCAsur4'};
            
        if strcmp(getenv('computername'), 'AVI-PC') 
            % laptop:
            allPruningFeatures = features_representative;
        else
            % workstation
            allPruningFeatures = features_smoothN;
            
        end
%         allFets = { ...
%          'Neg', 'PCAcuw3', 'PCAcuw4','PCAcuw5', 'PCAcuw6', 'PCAcuw7', 'PCAcuw8','PCAcuw10', 'PCAcuw12', 'PCAcuw16', ...
%                 'PCAcur3', 'PCAcur4','PCAcur5', 'PCAcur6', 'PCAcur7', 'PCAcur8','PCAcur10', 'PCAcur12', 'PCAcur16', ...
%                 'PCAsur1', 'PCAsur2', 'PCAsur3', 'PCAsur4'};
            
            
%         allFets = { ...
%         'GLFcuw4', 'GLFcuw8', 'GLFcuw12', 'GLFcuw16', ...
%         'Neg', ...
%          'PCAcur4' };    
%         allFets = { ...
%         'Neg', 'GLFcuw4', 'GLFcuw8', 'PCAcuw4', 'PCAcuw8'};    
    
     nPruningSets = length(allPruningFeatures);
     fprintf('Running "processDatafile" for the following %d pruning feature sets: \n', nPruningSets);
     fprintf(['   ' repmat('%s, ', 1, 9) '\n'], allPruningFeatures{:});
     fprintf('\n\n');

  
     CPU_id_arg = iff(isempty(CPU_id), 'all', num2str(CPU_id));
     for i = 1:nPruningSets
         fet = allPruningFeatures{i};
%          fet_str = allFet_str{i};
%          curSortingFeatures(fet);

         curPruningFeatures(fet);
         
         tmpfname = [spikePath sprintf('_workingOn_PruningSet_%s__%d_of_%d_(CPU_%s).tmp', fet, i, nPruningSets, CPU_id_arg)];
         save(tmpfname);         
         
         fprintf( '\n\n *** Pruning Feature set: %s  (%d/%d) ******* \n\n', fet, i, nPruningSets);
         
         processDatafile([], 'D'); 
         processDatafile([], 'N'); 
         
         delete(tmpfname);
         
         if exist([spikePath 'stop_running'], 'file')
             return;
         end
     end

     if strcmp(getenv('computername'), 'AVI-WORK-PC')
         s = dir([spikePath '_workingOn*']);
         if isempty(s)
             disp('Sending Email to self...');
             sendEmailToSelf('Completed all calculations');
         end
     end


end



%{
function processDatafile_helper

     if exist([spikePath 'stop_running'], 'file')
        movefile([spikePath 'stop_running'], [spikePath 'stop_runnin'])
     end
        
%     allSortingFeatures = {'Neg', 'PCAcuw8', 'PCAcuw6', 'PCAcuw4', 'Egy_PCAsnr1'};
%     allSortingFeatures = { {'Neg'}, {'PCAcuw6'}, {'Egy', 'PCAsnr1'}};
%     allSortingFeatures = { {} };

    features_list_S = dir('C:\ExperimentDB\Spikes\clusters');
    features_list = {features_list_S.name};
%     fet_strs = features_list( cellfun(@(s) isstrprop(s(1), 'alpha'), features_list ) );
    fet_strs = features_list( cellfun(@(s) strncmp(s, 'GLF', 3), features_list ) );
%     features_list = features_list(4:end-1);
    
%     allSortingFeatures = features_list; %fet_strs ;
    allSortingFeatures = fet_strs ;
    
    for i = 1:length(allSortingFeatures);
        fet = allSortingFeatures(i);
        curSortingFeatures(fet{:});
        
        fprintf( '\n\n *** Sorting Feature set: %s  ******* \n\n', fet{:});
        
        processDatafile;
        
        if exist([spikePath 'stop_running'], 'file')
            return;
        end           
    end








end
%}