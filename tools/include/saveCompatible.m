% saveCompatible.m

% Define library files that we want to standardise version on.
% Specify submodels of library files first to avoid warnings.
files = { ...
    'lib/cjt_library_controllers.mdl', ...
    'lib/cjt_library_controllers_experimental.mdl', ...
    'lib/cjt_library_models.mdl', ...
    'lib/cjt_library_observers.mdl', ...
    'lib/cjt_library.mdl', ...
    'examples/simulink/controlBlockDevelopment_generic.mdl', ...
    'examples/simulink/Ex_01_Open_Loop.mdl', ...
    'examples/simulink/Ex_02_Current_Control.mdl', ...
    'examples/simulink/Ex_03_I_Control_Plus_Feedforward.mdl', ...
    'examples/simulink/Ex_04_PD_Control_Plus_Feedforward.mdl', ...
    'examples/simulink/Ex_05_PD_Control.mdl', ...
    'examples/simulink/Ex_06_PD_Feedforward_OpenLoopDOB.mdl', ...
    'examples/simulink/Ex_07_PD_OpenLoopDOB.mdl', ...
    'examples/simulink/Ex_08_PD_Feedforward_ClosedLoopDOB.mdl', ...
    'examples/simulink/Ex_09_PD_ClosedLoopDOB.mdl', ...
    'examples/simulink/Ex_10_Open_Loop_With_Disturbance_Observer.mdl', ...
    'examples/simulink/Ex_11_Passivity_Based_Control.mdl', ...
    'examples/simulink/Ex_12_Passivity_Based_Control_Plus_Disturbance_Observer.mdl', ...
    'examples/simulink/openLoop.mdl', ...
    'examples/simulink/openLoop_inputFromFile.mdl', ...
    'unit_tests/simulinkBlockLibraryTest.mdl', ...
};

% Define Simulink versions (uncomment the one to save to)
%toSimulinkVersion       = '8.4';    % R2014b
%toSimulinkVersionName   = 'R2014b'; % R2014b
toSimulinkVersion       = '8.7';    % R2016a
toSimulinkVersionName   = 'R2016a'; % R2016a
%toSimulinkVersion       = '8.9';    % R2017a
%toSimulinkVersionName   = 'R2017a'; % R2017a
simulinkFileType        = 'MDL';    % Simulink file type (SLX or MDL)

% Counter for number of converted files
j = 0;

% Go through each file
for i=1:length(files)
    % Get filename
    file = files{i};
    
    % Check if file exists
    if (~exist(file, 'file'))
        error(['ERR: Library file ''' file ''' does not exist!']);
    end
    
    % Get file info and check SimulinkVersion (ReleaseName shows "<unknown
    % release>" for newer versions of MATLAB than the used one..)
    fileInfo = Simulink.MDLInfo(file);
    
    if (~strcmpi(fileInfo.SimulinkVersion, toSimulinkVersion))
        disp(['Saving library file ''' file ''' to Simulink ' toSimulinkVersion ' (MATLAB ' toSimulinkVersionName ')...']);
        
        % Get new filename/path
        [pathstr,name,ext] = fileparts(file);
        if (strcmp(pathstr, ''))
            newFile = [name '_compat' ext];
        else
            newFile = [pathstr filesep name '_compat' ext];
        end
        
        % Load old file and save standardised version to new filename
        load_system(file);
        Simulink.exportToVersion(name, newFile, [toSimulinkVersionName '_' simulinkFileType]);
        close_system(file);
        
        % Delete old file and move new file into its place
        delete(file);
        status = movefile(newFile, file);
        if (status ~= 1)
            error(['ERR: Could not move file ''' newFile ''' to ''' file '''']);
        end
        
        % Increment counter
        j = j + 1;
        
    else
        disp(['Library file ''' file ''' is OK.']);
    end
    
end

disp(['Finished. Converted ' num2str(j) ' files.']);

