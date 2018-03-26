% saveCompatible.m

% Define library files that we want to standardise version on.
% Specify submodels of library files first to avoid warnings.
files = { ...
    'lib/cjt_library_controllers.mdl', ...
    'lib/cjt_library_controllers_experimental.mdl', ...
    'lib/cjt_library_models.mdl', ...
    'lib/cjt_library_observers.mdl', ...
    'lib/cjt_library.mdl', ...
    'controlled.mdl', ...
};

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
    
    if (~strcmpi(fileInfo.SimulinkVersion, '8.7'))
        disp(['Saving library file ''' file ''' to Simulink 8.7 (MATLAB R2016a)...']);
        
        % Get new filename/path
        [pathstr,name,ext] = fileparts(file);
        if (strcmp(pathstr, ''))
            newFile = [name '_2' ext];
        else
            newFile = [pathstr filesep name '_2' ext];
        end
        
        % Load old file and save standardised version to new filename
        load_system(file);
        Simulink.exportToVersion(name, newFile, 'R2016a_MDL');
        close_system(file);
        
        % Delete old file and move new file into its place
        delete(file);
        status = movefile(newFile, file);
        if (status ~= 1)
            error(['ERR: Could not move file ''' newFile ''' to ''' file '''']);
        end
    end
    
end

disp('Finished.');
