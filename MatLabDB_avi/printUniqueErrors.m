
allErrors = {};
vars = who('celldata*');
for vi = 1:length(vars)
    var_name = vars{vi};
    vr = eval(var_name);
    if ischar(vr)
        disp([var_name ': ' vr]);
        allErrors = {allErrors{:}, vr};
    elseif isstruct(vr) && isfield(vr, 'id')
        disp([var_name ': ' vr.id ' : ' vr.txt]);
        allErrors = {allErrors{:}, vr.id};
    end
   
%     [uniqueErrors indices] = uniqueCount(allErrors);
%     ns = cellfun(@length, indices, 'UniformOutput', false);
%     {b(:); ns{:} }'
end
format long
% a(:,1) = cellfun(@length, indices, 'UniformOutput', false)';
% a(:,2) = uniqueErrors';
varBreakdown(allErrors);



% for MOVIES:
%     [  2]    'db:fps'                
%     [ 23]    'db:invalidTbTe'        
%     [  8]    'db:syncNframesMismatch'
%     [ 16]    'db:syncsTbTeMismatch'  
%     [217]    'file:invalidMovieFile'

% for NOISE
%     [  2]    'db:invalidTbTe'        
%     [ 31]    'db:syncNframesMismatch'
%     [ 60]    'db:syncsTbTeMismatch'  
%     [285]    'stim:interrupted'      

