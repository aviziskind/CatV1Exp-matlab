function printCellAcceptanceCriteria(cmpType, subtractSpont_flag, timeWindow)

%     printDegreeTuningCellsStats('all');

    if (nargin < 1) || isempty(cmpType)
        cmpType = curCmpType('');
    end

    if (nargin < 2) || isempty(subtractSpont_flag);
        subtractSpont = curSubtractSpont;
    else
        subtractSpont = exist('subtractSpont_flag', 'var') && isequal(subtractSpont_flag, 1);
    end
    subtractSpont_str = iff(subtractSpont, '(Spontaneous subtracted)', '(Spontaneous included)');
    
    useSigBckgResponseForSS = 0;
    
    if (nargin < 3) || isempty(timeWindow);
        timeWindow = curTimeWindow;
    end

    criteria_fn_d = getFileName('cellSelectedStats', [], [], struct('gratingType', 'drifting', 'cmpType', cmpType, 'subtractSpont', subtractSpont, 'timeWindow', timeWindow));
    criteria_fn_f = getFileName('cellSelectedStats', [], [], struct('gratingType', 'flashed', 'cmpType', cmpType, 'subtractSpont', subtractSpont, 'timeWindow', timeWindow));
    S_f = load(criteria_fn_f);
    S_d = load(criteria_fn_d);
    

        switch cmpType

            case 'degree',
               
%                 if ~subtractSpont
                    S_ori_f = S_f.S_criteria_ori;
                    S_ori_d = S_d.S_criteria_ori;
                    
                    fprintf('\n*** ORIENTATION %s -- CELL SELECTION STATS ***\n', subtractSpont_str);
                    fprintf('From a total of %d (%d) cells responding to drifting (flashed) gratings,\n', S_ori_d.n_Total_ori, S_ori_f.n_Total_ori); 
                    fprintf('%d (%d) cells passed the orientation selectivity criterion, \n', S_ori_d.n_Selective_ori, S_ori_f.n_Selective_ori);
                    fprintf('%d (%d) cells passed the reproducibility criterion, and \n', S_ori_d.n_Reproducible_ori, S_ori_f.n_Reproducible_ori);
                    fprintf('%d (%d) passed both of these criteria. \n', S_ori_d.n_Sel_Rep_ori, S_ori_f.n_Sel_Rep_ori);
                    fprintf('The orientation tuning curves of %d (%d) of these ... also had good fits ....\n', S_ori_d.n_Sel_Rep_Fit_ori, S_ori_f.n_Sel_Rep_Fit_ori);
                    pct_usedOri_d = S_ori_d.n_Used_ori / S_ori_d.n_Total_ori * 100;
                    pct_usedOri_f = S_ori_f.n_Used_ori / S_ori_f.n_Total_ori * 100;
                    fprintf('Thus, a total of %.1f%% (%.1f%%) of our cells passed all three selection criteria and were accepted for further analysis. \n\n', pct_usedOri_d, pct_usedOri_f);

                    S_spf_f = S_f.S_criteria_spf;
                    S_spf_d = S_d.S_criteria_spf;                    

                    fprintf('\n*** SPATIAL FREQUENCY %s -- CELL SELECTION STATS ***\n ', subtractSpont_str);

                    fprintf('From the %d (%d) cells whose spatial frequency tuning in response to drifting (flashed) gratings was studied, \n', S_spf_d.n_total_spf, S_spf_f.n_total_spf);
                    fprintf('%d (%d) had tuning curves that were reproducible, and\n', S_spf_d.n_reproducible_spf, S_spf_f.n_reproducible_spf);
                    fprintf('%d (%d) of these reproducible tuning curves were well fit by a skewed lognormal function. Thus, a total of \n', S_spf_d.n_rep_fit_spf, S_spf_f.n_rep_fit_spf);
                    pct_usedSpf_d = S_ori_d.n_Used_ori / S_ori_d.n_Total_ori * 100;
                    pct_usedSpf_f = S_ori_f.n_Used_ori / S_ori_f.n_Total_ori * 100;
                    fprintf('%.1f%% (%.1f%%) of cells recorded passed our spatial frequency tuning criteria.\n\n', pct_usedSpf_d, pct_usedSpf_f);

%                 else
%                     S_ori_d = S_d.S_criteria_ori;
%                     S_ori_f = S_f.S_criteria_ori;                    
%                     
%                     fprintf('\n*** ORIENTATION TUNING -- SPONTANEOUS SUBTRACTED ***\n');
%                     
%                     fprintf('From a total of %d (%d) cells responding to drifting (flashed) gratings,\n', S_ori_d.nOriTotal, S_ori_f.nOriTotal ); 
%                     fprintf('%d (%d) cells, or %.1f% (%.1f%), had responses significantly higher than spontaneous levels  according to this test.\n', ...
%                         S_ori_d.nOriSigResponse, S_ori_f.nOriSigResponse, S_ori_d.nOriSigResponse/S_ori_d.nOriTotal*100, S_ori_f.nOriSigResponse/S_ori_f.nOriTotal*100);
%                     fprintf('When combined with the original three acceptance criteria (significant orientation tuning, \n');
%                     fprintf('703 (287), or 64.4% (53.4%) satisfied all four criteria. This is compared to the \n', S_ori_d.nOriSel_Rep_Fit_Sig, S_ori_f.nOriSel_Rep_Fit_Sig, S_ori_d.nOriSel_Rep_Fit_Sig/S_ori_d.nOriTotal*100, S_ori_f.nOriSel_Rep_Fit_Sig/S_ori_d.nOriTotal*100);
%                     fprintf('761 (337) cells, or 69.8% (62.8%) that satisfied the original three criteria.\n', S_ori_d.nOriSigResponse, S_ori_f.nOriSigResponse);
%                     
%                 end
                
            case 'phase',
                fprintf('*** PHASE TUNING CURVES -- CELL SELECTION STATS *** ');

                fprintf('A total of 589/650 (90.6%) cells recorded during the drifting gratings and\n');
                fprintf('419/537 (78%) of the cells recorded during the flashed grating experiments passed this reproducibility test.\n\n');

                fprintf('Of the 419 cells with reproducible responses,\n');
                fprintf('149 (36%) had reproducible MIDs. This low percentage is largely due to the complex cells, ... :\n');
                fprintf('only 55/266 (21%) passed our criterion. In comparison,\n');
                fprintf('82% (94/114) of the MIDs from simple cells were reproducible. \n');
                fprintf('... We accepted the Gabor fit to the MID if the R2 goodness of fit measure was above our threshold of 0.5. \n');
                fprintf('A total of 109 (73%) of the Gabor fits satisfied this criterion.\n');


        end
end