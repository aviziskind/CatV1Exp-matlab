function t = durationOfExp(s)

    st = getGratingStimType( s.Gid );
    
    t = s.tempPeriod_sec * (st.nTrials * st.nOri * st.nSpf);
    3;

end