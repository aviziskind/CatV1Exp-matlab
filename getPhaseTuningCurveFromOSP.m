function y = getPhaseTuningCurveFromOSP(OSP)

    if isstruct(OSP)
        [OSP, oris, spfreqs, phases] = elements(OSP);
    end

    OS = mean(OSP, 3);
    [tmp, inds] = maxElement(OS);
    [ori_ind, sp_ind] = elements(inds);

    y = squeeze(OSP(ori_ind, sp_ind, :));
end