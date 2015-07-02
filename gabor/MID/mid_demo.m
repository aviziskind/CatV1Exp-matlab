function mid_demo
    
    cellName = 'democell';
    dataDir = [baseDir 'Demo\My_Data'];
    stimFileName = [baseDir 'Demo\walk5_ieee_128x128x16384.raw'];
    spikeFileName = [baseDir 'Demo\demo_spike_Noise0.50_Rep1_1.isk'];
    dimx = 128; dimy = 128;
    nFrames = 16384;
    x0 = 5; y0 = 5;
    dx = 120; dy = 120;
    cx = 4; cy = 4;
    nLags = 1;                
    funcName = 'do_MaxInfoDim_Nvec';
    
%         do_MaxInfoDim_Nvec democell
%         Demo/walk5_ieee_128x128x16384.raw 
%         Demo/demo_spike_Noise0.50_Rep1_1.isk 
%         Demo/My_Data 128 128 16384 5 5 120 120 4 4 1 32 2
                
    cmd_str = sprintf('%s %s %s %s %s %d %d %d %d %d %d %d %d %d %d %d', [baseDir funcName], cellName, ...
        stimFileName, spikeFileName, dataDir, dimx, dimy, nFrames, x0, y0, dx, dy, cx, cy, nLags);

    dos(cmd_str);   
    
    dimInfo = [round(x0+dx-1)/cx, round(y0+dy-1)/cy, nLags, cx];
    MIDfilename  = [dataDir 'Test_v_' cellName];


    v_MID = read_vec_pxpxt_file (MIDfilename,dimInfo);
    
    MID = mean(v_MID{1}, 3);
    
    
    figure(10); imagesc(MID);
end

%     dimInfo = [dx/cx, dy/cx, nLags, cx];
%     Test_v_Group4470_Cell0_16x16x1_4_jack1.dat
%     Test_v_Group4470_Cell0_64x64x1_4_jack1.dat
