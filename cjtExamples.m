function cjtExamples


close all;

figname = 'Compliant Joint Toolbox Examples';

fig_props = { ...
    'name'                   figname ...
    'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            'normal' ...
    'visible'                'on' ...
    'units'                  'normalized' ...
    'createfcn'              ''    ...
    'closerequestfcn'        'delete(gcbf)' ...
    };

% 'CreateFcn'              @localCreate_callback

% Create figure to use as GUI in your main function or a local function
mainFig = figure(fig_props{:});
% create structure of handles
fig_handles = guihandles(mainFig);

cjtPath = fileparts(which('cjtExamples'));
examplePath = [cjtPath, filesep, 'examples'];


%% Gater information about what examples are available
matlabExamples = scanForExamples([examplePath, filesep, 'matlab']);
simulinkExamples = scanForExamples([examplePath, filesep, 'simulink']);

nMatlabExamples = numel(matlabExamples);
matlabExampleStrings = cell(nMatlabExamples,1);
for iEx = 1:nMatlabExamples
    matlabExampleStrings{iEx,1} = matlabExamples(iEx).displayName;
end

nSimulinkExamples = numel(simulinkExamples);
simulinkExampleStrings = cell(nSimulinkExamples,1);
for iEx = 1:nSimulinkExamples
    simulinkExampleStrings{iEx,1} = simulinkExamples(iEx).displayName;
end


%% Some parameters to position and align uicontrols
listHeight = 0.3;
listWidth = 0.4;
listX = 0.05;
listY = 0.6;
btnY = 0.03;
btnWidth = 0.25;

%% Uicontrols for Matlab Examples
fig_handles.matlabTitle = uicontrol('Style','text',...
    'units','normalized',...
    'String','Matlab Examples',...
    'FontWeight', 'bold',...
    'FontSize', 12,...
    'HorizontalAlignment', 'left',...
    'Tag','matlabTitle_txt',...
    'Position', [listX 0.91 listWidth 0.05],...
    'Callback', {@localListbox_callback});

fig_handles.matlabListbox = uicontrol('Style','listbox',...
    'units','normalized',...
    'String',matlabExampleStrings,...
    'BackgroundColor','w',...
    'Tag','matlabListbox',...
    'Position',[listX listY listWidth listHeight ],...
    'Callback', {@localListbox_callback});

%% Uicontrols Simulink Examples
fig_handles.simulinkTitle = uicontrol('Style','text',...
    'units','normalized',...
    'String','Simulink Examples',...
    'FontWeight', 'bold',...
    'FontSize', 12,...
    'HorizontalAlignment', 'left',...
    'Tag','simulinkTitle_txt',...
    'Position', [0.5+listX 0.91 listWidth 0.05],...
    'Callback', {@localListbox_callback});


fig_handles.simulinkListbox = uicontrol('Style','listbox',...
    'units','normalized',...
    'String',simulinkExampleStrings,...
    'BackgroundColor','w',...
    'Tag','simulinkListbox',...
    'Position',[0.5+listX listY listWidth listHeight ],...
    'Callback', {@localListbox_callback});

%% Common uicontrols
% Run an example
fig_handles.run_btn = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'String','Run Example',...
    'Tag','run_btn',...
    'Position',[listX btnY btnWidth 0.05],...
    'Callback',{@localRun_callback});

% Open the example m-file
fig_handles.open_btn = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'String','Open Example',...
    'Tag','open_btn',...
    'Position',[listX+btnWidth+0.05 btnY btnWidth 0.05],...
    'Callback',{@localOpen_callback});

% Open the example m-file
fig_handles.close_btn = uicontrol('Style','pushbutton',...
    'units','normalized',...
    'String','Close Dialog',...
    'Tag','close_btn',...
    'Position',[listX+2*(btnWidth+0.05) btnY btnWidth 0.05],...
    'Callback',{@localClose_callback});

descriptionString = '';
fig_handles.descriptionListbox = uicontrol('Style','listbox',...
    'units','normalized',...
    'String',descriptionString,...
    'Tag','matlabListbox',...
    'Position',[listX 0.1 0.9 0.43 ],...
    'Callback', {@localDescription_callback});

fig_handles.descriptionTitle = uicontrol('Style','text',...
    'units','normalized',...
    'String','Example Description',...
    'FontWeight', 'bold',...
    'FontSize', 12,...
    'HorizontalAlignment', 'left',...
    'Tag','description_txt',...
    'Position', [listX 0.54 listWidth 0.05]);%,...

% Extend the handles struct by some information that might be useful in
% callbacks.
fig_handles.cjtPath = cjtPath;
fig_handles.examplePath = examplePath;
fig_handles.matlabExamples = matlabExamples;
fig_handles.simulinkExamples = simulinkExamples;
fig_handles.mainFig = mainFig;

% Save the structure
guidata(mainFig,fig_handles);

end

function localDescription_callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns contents
% contents{get(hObject,'Value')} returns selected item from listbox1

% Get the structure using guidata in the local function
fig_handles = guidata(gcbo);

% Find out which example has been selected.
listboxTag = get(fig_handles.activeListbox,'Tag');
index_selected = get(fig_handles.activeListbox,'Value');

% Get m filename
switch listboxTag
    case 'matlabListbox'
        mFileName = fig_handles.matlabExamples(index_selected).fileName;
    case 'simulinkListbox'
        mFileName = fig_handles.simulinkExamples(index_selected).fileName;
end

% Extract example description
descriptionText = extractExampleDescription(mFileName);

set(fig_handles.descriptionListbox,'string', descriptionText);

end


function localClose_callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns contents
% contents{get(hObject,'Value')} returns selected item from listbox1

% Get the structure using guidata in the local function
fig_handles = guidata(gcbo);

close(fig_handles.mainFig);


end

function localOpen_callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns contents
% contents{get(hObject,'Value')} returns selected item from listbox1

% Get the structure using guidata in the local function
fig_handles = guidata(gcbo);

% Find out which example has been selected.
listboxTag = get(fig_handles.activeListbox,'Tag');
index_selected = get(fig_handles.activeListbox,'Value');

% Get m filename
switch listboxTag
    case 'matlabListbox'
        mFileName = fig_handles.matlabExamples(index_selected).fileName;
    case 'simulinkListbox'
        mFileName = fig_handles.simulinkExamples(index_selected).fileName;
end

% Open in editor
open(mFileName);

end


function localListbox_callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns contents
% contents{get(hObject,'Value')} returns selected item from listbox1


% Get the structure using guidata in the local function
fig_handles = guidata(gcbo);

% Update the active listbox property.
fig_handles.activeListbox = hObject;

% If a double-click has occurred, directly run the example.
if strcmp(get(fig_handles.mainFig,'SelectionType'),'open')
    disp('Double Click detected!')
    localRun_callback(hObject, [], [])
end

guidata(hObject,fig_handles);

% Update example description field.
localDescription_callback(hObject, [], []);

end



function [examples] = scanForExamples(examplePath)

examples = struct([]);

% First scan the matlab examples
dirContents = dir([examplePath,filesep,'*.m']);

nFiles = numel(dirContents);

for iFile = 1:nFiles
    
    fName = [examplePath, filesep, dirContents(iFile).name];
    
    dispName = extractExampleDisplayName(fName);
    
    if ~isempty(dispName)
        localExample.fileName = fName;
        localExample.displayName = dispName;
        examples = [examples, localExample];
    end
    
end


end

function dispName = extractExampleDisplayName(fName)

dispName = [];

fileContents = fileread(fName);

[startIndex,endIndex] = regexp(fileContents,'% #! [^\f\n\r\t\v]*[\f\n\r\t\v]');

dispName = fileContents(startIndex+5:endIndex);

end

function description = extractExampleDescription(fName)

description = {};

fid = fopen(fName);

stopFlag = 0;
cnt = 1;

while(~stopFlag)
    tline = fgetl(fid);
    
    if ~isempty(tline) && tline(1) == '%'
        description{cnt,1} = tline(3:end);
    else
        stopFlag = 1;
    end
    
    cnt = cnt +1;
    
end

fclose(fid);

end
