% S1 = load('psthWindowData_st_mean_allStats');
% S2 = load('psthWindowData_st_mean');
% 
% fn1 = fieldnames(S1);
% for i = 1:length(fn1)   
%     v1 = S1.(fn1{i});
%     v2 = S2.(fn1{i});
%     
%     v1.rep_p = v2.rep_p;
%     S1.(fn1{i}) = v1;
% end
% save('psthWindowData_st_mean_allStats', '-struct', 'S1', '-v6');

%%%%%%%%%%%%%%%%%

% S2 = load('psthWindowData_st_mean');

fn = fieldnames(S2);
for i = 1:length(fn)   
    v = S2.(fn{i});
    v = rmfield(v, {'tau_t'});
    S2.(fn{i}) = v;
end

% save('psthWindowData_st_mean', '-struct', 'S2', '-v6');