function [oris_shifted, otc_shifted, idx_shift] = shiftOriTuningCurve(oris, otc, ori_pref, ori_start)
    dOri = diff(oris(1:2));

    if nargin < 3
        ori_start = -90;
    end        
    
    oris_shifted = [0:dOri:360-dOri]+ori_start;
    otc_shft_amount = - ((oris_shifted(1)-oris(1))+ori_pref);
    
    nShift_actual = otc_shft_amount/dOri;
    nShift_round = round(nShift_actual);
    ori_shift_remain = -(nShift_round-nShift_actual)*dOri;
    
    oris_shifted = oris_shifted+ori_shift_remain;
    
    idx_shift = getIdxShift(nShift_round, length(oris));    
    otc_shifted = otc(idx_shift);
    
    3;
%     otc_shifted = 

end


function idx = getIdxShift(nShift, N)

    if nShift >= 0
        idx = [N-nShift+1:N, 1:N-nShift];
        
    elseif nShift < 0        
        idx = [-nShift+1:N, 1:-nShift];
        
    end
    assert(length(idx) == N);

end