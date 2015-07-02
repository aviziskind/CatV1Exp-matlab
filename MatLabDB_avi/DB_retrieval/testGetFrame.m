function testGetFrame
%     getMovieStimulusFrame('load', 'Gid', 2249)
%     frm = getMovieStimulusFrame(2053);

%     getNoiseStimulusFrame('load', 'presId', 60)
%     frm = getNoiseStimulusFrame(10);
    
    nframes = getGratingStimulusFrame('load', 'presId', 1479)
    for i = 1:1
        frm = getGratingStimulusFrame(i);
        imagesc(frm); colormap('gray'); colorbar; axis equal square tight xy; 
        fprintf('*');
        drawnow;
    end

%     getGratingStimulusFrame('load', 'Gid', 855)
%     frm = getGratingStimulusFrame(10);
%     
% f    3;
end

%     getGratingStimulusFrame('load', 'Did', 74)
%     for i = 1:250
%         frm = getGratingStimulusFrame(i);
%         imagesc(frm);
%         pause(1/120);
%         fprintf('*');
%         drawnow;
%     end
