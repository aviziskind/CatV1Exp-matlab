function [pair_passedAllCrit, general_str, specific_str] = getPairIdxsThatPassCriteria(criteria, measure_name, loc_i, pt_idxs, pairData, allMeasures)
    pair_passedAllCrit = true(length(pt_idxs), 1);
%%    
    pairData_fields = fieldnames(pairData);
    criteria_fields = fieldnames(criteria);
    general_crit_field_names = intersect(criteria_fields, pairData_fields);
    specific_criteria_field_names = intersect(criteria_fields, allMeasures);
    
    remaining_fields = setdiff(criteria_fields, [general_crit_field_names; specific_criteria_field_names]);
    assert(isempty(remaining_fields));
    
    general_criteria = struct;
    for fld_i = 1:length(general_crit_field_names)
        general_criteria.(general_crit_field_names{fld_i}) = criteria.(general_crit_field_names{fld_i});
    end
  %%      
    if isfield(criteria, measure_name)        
        specific_criteria = criteria.(measure_name);
    else 
        specific_criteria = struct;
    end
    gen_spec_criteria = {general_criteria, specific_criteria};
    nCriteriaTypes = length(gen_spec_criteria);
    criteria_strs = cell(1, nCriteriaTypes);
    
    
        %%
    for crit_type_i = 1:nCriteriaTypes  % 1: general (for all measures) 2: specific (only for this measure)
        criteria_type_i = gen_spec_criteria{crit_type_i};
        crit_flds = fieldnames(criteria_type_i);
        criteria_str = '';
        for crit_i = 1:length(crit_flds)
            crit_str = crit_flds{crit_i};
            crit = criteria_type_i.(crit_flds{crit_i});
            
            if isnumeric(crit)
                crit_value = crit;
                crit_op = @gt;
                
            elseif isstruct(crit)
                crit_value = crit.value;
                crit_op = crit.op;
            end
            
            criteria_str_i = getPairCriteriaStr(crit_flds{crit_i}, crit);
            criteria_str = appendToStr(criteria_str, criteria_str_i, '; '); 
            
            pairData_vals = pairData.(crit_str)(pt_idxs,:,:);
            if strncmp(crit_str, 'loc_', 4);
                pairData_vals = pairData_vals(:, loc_i); % choose data for this specific location
            end    
            
            passedThisCrit = crit_op(pairData_vals, crit_value);
            if ~any(passedThisCrit)
                3;
            end
            pair_passedAllCrit = pair_passedAllCrit & passedThisCrit;
        end
        criteria_strs{crit_type_i} = criteria_str;
    end
    
    general_str = criteria_strs{1};
    specific_str = criteria_strs{2};
end

