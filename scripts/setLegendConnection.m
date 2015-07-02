function setLegendConnection(h, newstatus)
    % newstatus = 'on' or 'off';
    set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle', newstatus);   
end