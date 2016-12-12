function cjtExamples

close all;

optStruct = struct;
optStruct.cjtPath = fileparts(which('cjtExamples'));
optStruct.examplePath = [optStruct.cjtPath, filesep, 'examples'];


%% Gater information about what examples are available
[matlabExamples, simulinkExamples] = scanForExamples(optStruct);

% hMainFigure = figure(...       % The main GUI figure
%                     'Name','Compliant Joint Toolbox Examples',...
%                     'MenuBar','none', ...
%                     'Toolbar','none', ...
%                     'HandleVisibility','callback', ...
%                     'Color', get(0,...
%                              'defaultuicontrolbackgroundcolor'));
                         
figname = 'Compliant Joint Toolbox Examples';

% listsize = (numel(matlabExamples))
% 
% fus = 8;
% ffs = 8;
% uh = 22;
% 
% promptstring = ['1'];
% 
% ex = get(0,'DefaultUicontrolFontSize')*1.7;  % height extent per line of uicontrol text (approx)
% 
% fp = get(0,'DefaultFigurePosition');
% w = 2*(fus+ffs)+listsize(1);
% h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh+(smode==2)*(fus+uh);
% fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed


fig_props = { ...
    'name'                   figname ...
    'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            'modal' ...
    'visible'                'on' ...
    'createfcn'              ''    ...
    'closerequestfcn'        'delete(gcbf)' ...
            };
        
liststring = '';

mainFig = figure(fig_props{:});

listbox = uicontrol('Style','listbox',...
                    'String',liststring,...
                    'BackgroundColor','w',...
                    'Tag','listbox',...
                    'Callback', {@doListboxClick});

% ok_btn = uicontrol('Style','pushbutton',...
%                    'String',okstring,...
%                    'Tag','ok_btn',...
%                    'Callback',{@doOK,listbox});
% 
% cancel_btn = uicontrol('Style','pushbutton',...
%                        'String',cancelstring,...
%                        'Tag','cancel_btn',...
%                        'Callback',{@doCancel,listbox});

% listbox = uicontrol('Style','listbox',...
%                     'Position',[ffs+fus ffs+uh+4*fus+(smode==2)*(fus+uh) listsize],...
%                     'String',liststring,...
%                     'BackgroundColor','w',...
%                     'Max',smode,...
%                     'Tag','listbox',...
%                     'Value',initialvalue, ...
%                     'Callback', {@doListboxClick});
% 
% ok_btn = uicontrol('Style','pushbutton',...
%                    'String',okstring,...
%                    'Position',[ffs+fus ffs+fus btn_wid uh],...
%                    'Tag','ok_btn',...
%                    'Callback',{@doOK,listbox});
% 
% cancel_btn = uicontrol('Style','pushbutton',...
%                        'String',cancelstring,...
%                        'Position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
%                        'Tag','cancel_btn',...
%                        'Callback',{@doCancel,listbox});

%     gui_Singleton = 1;
%     gui_State = struct(...
%         'gui_Name',       'Compliant Joint Toolbox Examples', ...
%         'gui_Singleton',  1, ...
%         'gui_OpeningFcn', @local_gui_OpeningFcn, ...
%         'gui_OutputFcn',  @local_gui_OutputFcn, ...
%         'gui_LayoutFcn',  [], ...
%         'gui_Callback',   []);
%     h = gui_mainfcn(gui_State);


end


% --- Executes just before rtbdemo_gui is made visible.
function local_gui_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to rtbdemo_gui (see VARARGIN)
    
    % Choose default command line output for rtbdemo_gui
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    
    local_initialize_gui(hObject, handles, false);
end

% --- Outputs from this function are returned to the command line.
function varargout = local_gui_OutputFcn(hObject, eventdata, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

function local_initialize_gui(fig_handle, handles, isreset)
    % Update handles structure
    guidata(handles.figure1, handles);
end

function [matlabExamples, simulinkExamples] = scanForExamples(optStruct)
    
matlabExamples = struct([]);
simulinkExamples = struct([]);

% First scan the matlab examples
mPath = [optStruct.examplePath, filesep ,'matlab'];
dirContents = dir([mPath,filesep,'*.m']);

nFiles = numel(dirContents);

for iFile = 1:nFiles

    fName = [mPath, filesep, dirContents(iFile).name];
    
    dispName = extractExampleDisplayName(fName);
    
    if ~isempty(dispName)
        localExample.fileName = fName;
        localExample.displayName = dispName;
        matlabExamples = [matlabExamples, localExample];
    end
    
end

% Then scan the simulink examples
sPath = [optStruct.examplePath, filesep ,'simulink'];
dirContents = dir([sPath,filesep,'*.m']);

    
end

function dispName = extractExampleDisplayName(fName)
    
    dispName = [];
    
    fileContents = fileread(fName);
    
    [startIndex,endIndex] = regexp(fileContents,'% #! [^\f\n\r\t\v]*[\f\n\r\t\v]');
    
    dispName = fileContents(startIndex+5:endIndex);

end