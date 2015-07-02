% curTimeWindow

responseTypes = [1, 2];
gratingTypes = [1, 2];
% phaseOE_modes = {'aa', 'oe_diff_k12'};
phaseOE_modes = {'oe_diff_k12'};

for ri = 1:length(responseTypes)
    curResponseType(responseTypes(ri));
    
    for gi = 1:length(gratingTypes)
        curGratingType(gratingTypes(gi));    
       
        for p_i = 1:length(phaseOE_modes);
            curPhaseOEmode(phaseOE_modes{p_i});
            gcmp;
%             if p_i == 1
%                 gen3;
%             elseif p_i > 1
%                 gcmp;
%             end
        end
    end
    
end

% curResponseType(1); curGratingType(1); gen3; curGratingType(2); gen3; curResponseType(2); curGratingType(1); gen3; curGratingType(2); gen3;
% curPhaseOEmode('oe_diff_k12')