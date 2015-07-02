function [hyperplaneCriteron_id_out, hyperplaneCriteron_out] = curHyperPlaneCriterion(s)

    persistent hyperplaneCriteron

    filename = [CatV1Path 'MatLabDB_avi' filesep 'curHPcriteronFile.mat'];
    
    setType = (nargin == 1) && ~isempty(s);
    
    if setType    % set current grating type 
    
        if isequal(s, 1) || strncmp(s, 'D', 1)
            hyperplaneCriteron = 'D';
        elseif isequal(s, 2) || strncmp(s, 'N', 1)
            hyperplaneCriteron = 'N';
        else
            error('invalid grating type: enter "D" or "N"');
        end
        hyperplaneCriteron_out = hyperplaneCriteron;
        save(filename, 'hyperplaneCriteron');
        
    else
        if isempty(hyperplaneCriteron)            
            if exist(filename, 'file')
                S = load(filename);
                hyperplaneCriteron = S.hyperplaneCriteron;
            else
                error('Hyperplane Criteron file does not exist');
            end
        end
        hyperplaneCriteron_out = hyperplaneCriteron;

        allhyperplaneCriterons = {'D', 'N'};
        hyperplaneCriteron_id_out = find(strcmp(hyperplaneCriteron, allhyperplaneCriterons)); 
        
        if (nargin==1) && isempty(s)
            hyperplaneCriteron_id_out = hyperplaneCriteron_out;
        end                
    end

end
