function C = adjustVariableDialog(varId, elId)

        curMin = vMin{varId}(elId);
        curMax = vMax{varId}(elId);
        curStepSize = vStep{varId}(elId);
        curVal = varValues{varId}(elId);
        curNsteps = (curMax-curMin)/curStepSize;
        varName = varNames{varId};


hDlg = dialog('Name','adjustVariableRange',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[520 328 341 172],...
'Resize','off',...
'WindowStyle',get(0,'defaultfigureWindowStyle'),...
'Visible','on');

hOK     = uicontrol('Parent',hDlg, 'Position',[250 108 69 30], 'String','OK' );
hCancel = uicontrol('Parent',hDlg, 'Position',[250 71 69 30], 'String','Cancel' );

uicontrol('Style','text', 'String','Adjust Range for variable:', 'Parent',hDlg, 'FontSize',11, 'Position',[78 150 190 18]);
uicontrol('Style','text', 'String','Max', 'Parent',hDlg, 'HorizontalAlignment','right', 'Position',[0 83 70 20] );
uicontrol('Style','text', 'String','Min', 'Parent',hDlg, 'HorizontalAlignment','right', 'Position',[0 115 70 20]);
uicontrol('Style','text', 'String','Step size', 'Parent',hDlg, 'HorizontalAlignment','right', 'Position',[0 51 70 20]);
uicontrol('Style','text', 'String','Num Steps', 'Parent',hDlg, 'HorizontalAlignment','right', 'Position',[0 17 70 20]);

curMin

hMinField = uicontrol('Style','edit', 'String', num2str(vMin{varId}(elId)), 'Parent',hDlg, 'Position',[85 119 100 20], 'BackgroundColor',[1 1 1]);
hMaxField = uicontrol('Style','edit', 'String', num2str(vMax{varId}(elId)), 'Parent',hDlg, 'Position',[85 87 100 20], 'BackgroundColor',[1 1 1], 'Tag','edit3');
hNumStepsField = uicontrol('Style','edit', 'String','dNumSteps', 'Parent',hDlg, 'Position',[85 21 100 20], 'BackgroundColor',[1 1 1], 'Tag','edit4');
hStepSizeField = uicontrol('Style','edit', 'String', num2str(vMin{varId}(elId)), 'Parent',hDlg, 'Position',[85 55 100 20], 'BackgroundColor',[1 1 1], 'Tag','edit5');
hSlider = uicontrol('Style','slider', 'Parent',hDlg, 'Position',[186 21 18 19], 'BackgroundColor',[1 1 1], 'Value', 10, 'Min', 1, 'Max', 1000);



end