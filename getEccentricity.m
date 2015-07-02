function e = getEccentricity(Gid)

    fovea_offset_deg = [17;7];

    sd = siteDataFor(Gid);
    sd.locationData.LocId;

    % get Optic Disk location
    [leftDiscLoc_pix, rightDiscLoc_pix, locationHemi] = getOpticDiskLocation(idType, idVal, leftRight);
    degPerPix = sd.stimulusInfo.degreesPerBlock/sd.stimulusInfo.pixPerBlock;
        
    leftDiscLoc_deg  = leftDiscLoc_pix*degPerPix;
    rightDiscLoc_deg = rightDiscLoc_pix*degPerPix;
    
    leftHemi = strcmp(locationHemi, 'L');
    rightHemi = strcmp(locationHemi, 'R');
    
    useLeftDisc = (leftHemi && ipsiEye) || (rightHemi && contraEye);
    useRightDisc = (rightHemi && ipsiEye) || (leftHemi && contraEye);
    
    if useLeftDisc
        fovea_deg = [leftDiscLoc_deg + fovea_offset_deg(1);
                     leftDiscLoc_deg + fovea_offset_deg(2)];
        
    elseif useRightDisc
        fovea_deg = [rightDiscLoc_deg-fovea_offset_deg(1);
                     rightDiscLoc_deg+fovea_offset_deg(2)];
        
    end
    
    stimulusCenter_deg = dbGetStimulusCenter(Gid);
    
    ecc = norm( fovea_deg - stimulusCenter_deg );
    
end







% function calcGroupEccentricity(Gid)
% 
% 
% 
% 
% 
% 
% end