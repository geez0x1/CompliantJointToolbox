classdef jointBuilder
    %JOINTBUILDER Builds joint model files
    
    properties
        buildDir = ['.',filesep,'build'];
        overwrite = 0;
    end
    
    methods
        %__________________________________________________________________
        % Constructor
        function this = jointBuilder()
            
        end
        
        %__________________________________________________________________
        % Remove all created files
        function purge(this)
            rmdir(this.buildDir,'s');
        end
        
        %__________________________________________________________________
        % Build a joint object file
        function className = buildJoint(this, paramName, modelName, nonlinearModelName)
            % buildJoint(paramName, modelName, nonlinearModelName)
            % Builds a joint model file in the build/ directory based on
            % the supplied parameter and model names
            
            % Process nonlinear model names
            if ~exist('nonlinearModelName', 'var')
                nonlinearModelName = '';
                nonlinearModelNameCellStr = '''''';
            else
                nonlinearModelName_multiple = 0; % default
                if (iscell(nonlinearModelName) && length(nonlinearModelName) > 1)
                    nonlinearModelName_multiple = 1;
                    nonlinearModelNameStr = strjoin(nonlinearModelName, '_');
                    nonlinearModelNameCellStr = strcat('{''', strtrim(strjoin(nonlinearModelName, ''', ''')), '''}');
                else
                    nonlinearModelNameStr = nonlinearModelName;
                    nonlinearModelNameCellStr = ['''' nonlinearModelName ''''];
                end
            end
            
            % Create class name
            eval(paramName);
            className = [paramName '_' modelName];
            if ~isempty(nonlinearModelName)
                className = [className '_' nonlinearModelNameStr];
            end
            className(isspace(className)) = [];
            
            % Fix long filenames
            classNameLong = className;
            if (length(className) > 63-2)
                className = [className(1:63-2-3) '_xx'];
            end
            
            % Joint name is equal to the long class name
            jointName = classNameLong;
            
            % Create build directory if necessary
            if ~exist(this.buildDir,'dir')
                mkdir(this.buildDir);
            end
            
            % Create class filename
            classFName = [this.buildDir,filesep,className,'.m'];
            
            % Confirm file overwrite
            fid = [];
            button = 'yes';
            if exist(classFName,'file') && this.overwrite == 0
                button = questdlg(sprintf('%s exists! Overwrite?',classFName),'Overwrite file?');               
            end
            switch lower(button)
                case {'no', 'cancel'}
                    fprintf('File %s exists. Aborted.\n', classFName);
                    return;
                case 'yes'
                    fid = fopen(classFName,'w+');
            end
            
            %% Class definition
            fprintf(fid,'%% %s\n\n', classNameLong);
            fprintf(fid,'classdef %s < genericJoint\n\n', className);
            
            % Description
            % ...
            
            %% Private properties
            fprintf(fid,'\tproperties (SetAccess = private)\n');
            fprintf(fid,'\tend\n\n');
            
            %% Methods
            fprintf(fid,'\tmethods\n');
            
            % Constructor
            paramsStr = fileread([paramName '.m']);
            paramsStr = regexprep(paramsStr, '^', '\t', 'emptymatch', 'lineanchors');
            constructorStr = fileread('constructor.m');
            constructorStr = sprintf(constructorStr, className, paramsStr, jointName, paramName, modelName, nonlinearModelNameCellStr);
            constructorStr = regexprep(constructorStr, '^', '\t\t', 'emptymatch', 'lineanchors');
            constructorStr = regexprep(constructorStr, '%', '%%', 'emptymatch', 'lineanchors');
            fprintf(fid, [constructorStr '\n\n']);
            
            % getDynamicsMatrices
            getDynStr = fileread('getDynamicsMatrices.m');
            getDynStr = regexprep(getDynStr, '^', '\t\t', 'emptymatch', 'lineanchors');
            fprintf(fid, [getDynStr '\n\n'], modelName);
            
            % getNonlinearDynamics
            if (~isempty(nonlinearModelName) && nonlinearModelName_multiple)
                % Multiple nonlinear terms
                getNonlinDynStr = fileread('getNonlinearDynamics_multiple.m');
                fStr = '';
                for i=1:length(nonlinearModelName)
                   fStr = strcat(fStr, nonlinearModelName(i), '(obj, x, dx)');
                   if (i < length(nonlinearModelName))
                       fStr = strcat(fStr, '+');
                   end
                end
                getNonlinDynStr = regexprep(getNonlinDynStr, '^', '\t\t', 'emptymatch', 'lineanchors');
                fprintf(fid, [getNonlinDynStr '\n\n'], char(fStr));
            elseif ~isempty(nonlinearModelName)
                % Single nonlinear term
                getNonlinDynStr = fileread('getNonlinearDynamics.m');
                getNonlinDynStr = regexprep(getNonlinDynStr, '^', '\t\t', 'emptymatch', 'lineanchors');
                fprintf(fid, [getNonlinDynStr '\n\n'], nonlinearModelName);
            else
                % No nonlinear terms
                getNonlinDynStr = fileread('getNonlinearDynamics_none.m');
                getNonlinDynStr = regexprep(getNonlinDynStr, '^', '\t\t', 'emptymatch', 'lineanchors');
                fprintf(fid, [getNonlinDynStr '\n\n'], nonlinearModelName);
            end
            
            % End methods
            fprintf(fid,'\tend\n\n');
            
            %% End class definition
            fprintf(fid,'end\n');
            
            %% Close file
            fclose(fid);
            
            % Status
            disp(['Joint built and written to ' classFName]);
        end
        
    end
    
end

