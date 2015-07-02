%go through all values
global statVarHandles 
% allStatNames = {'cc_p', 'rho_p', 'rho_p_nz', 'rho_p_nznz', 'r_entropy'};
allStatNames = {'cc_p', 'rho_p', 'rho_p_nznz', 'r_entropy'};
intSizes = 3:9;
nStdThs = 2:.1:4;
objThs = 2:.1:4;

% st_i = 1; sz_i = 1;std_i = 1;

progressBar('init=', length(allStatNames)*length(intSizes)*length(nStdThs)*length(objThs))
for st_i = 1:length(allStatNames)
    for sz_i = 1:length(intSizes)
        for std_i = 1:length(nStdThs)
            for oi_i  = 1:length(objThs)
                fprintf('%d/%d, %d/%d, %d/%d, %d/%d ... ', st_i, length(allStatNames), sz_i, length(intSizes), std_i, length(nStdThs), oi_i, length(objThs) ); 
                tic;
                manipulateSet(statVarHandles, {'statName', 'intSize', 'nStdTh', 'objTh'}, ...
                    {allStatNames{st_i}, intSizes(sz_i), nStdThs(std_i), objThs(oi_i)});
                toc;
                progressBar;
            end
        end
    end
end

%              {'intSize', [1:20], intSize0}, ...             
%              {'nStdTh', [.1:.1:4], nStdTh0}, ...
%              {'objTh', [.1:.1:4], objTh0}, ...
