function criteria_str = getPairCriteriaStr(critName, critCondition)
        
    if isnumeric(critCondition)
        crit_value = critCondition;
        crit_op = @gt;

    elseif isstruct(critCondition)
        crit_value = critCondition.value;
        crit_op = critCondition.op;
    end
    op_str = op2str(crit_op);

    if strncmpi(critName, 'SCtype', 6)
        crit_value_str = sc_idx2str(crit_value);
    else
        crit_value_str = num2str(crit_value);
    end

    criteria_str = [critName ' ' op_str ' ' crit_value_str];

end
            

function sc_str = sc_idx2str(sc_idx)
    switch sc_idx
        case 0, sc_str  = 'C/C';
        case 1, sc_str  = 'S/C';
        case 2, sc_str  = 'S/S';
        otherwise , sc_str  = '?';
    end

end