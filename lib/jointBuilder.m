% JOINTBUILDER Class use to create joint models.
%
% Joint models consist of linear and nonlinear dynamics equations along
% with a parameter set. The number and complexity of the nonlinear effects 
% to be considered vary for a particular joint, depending on the
% respecrtive objective. The same principal dynamics lead to diverse
% results when parameterized differently. 
%
% This class allows to create dynamics models with arbitrary combinations
% of the provided dynamics templates and parameter sets. It automatically
% assembles those in a class file implementic a genericJoint.
%
% Notes::
%
%
% Examples::
%
%     %% Instantiate a jointBuilder
%     jb = jointBuilder;
% 
%     %% Build joint model classes
%     % A linear model of a WALKMAN leg 
%     % actuator with 10 kNm/rad torsion 
%     % bar stiffness 
%     jb.buildJoint('WMBig10k',... parameters
%         [],... full linear dynamics
%         [],... no nonlinear dynamics
%         [],... static electric subsystem
%         'my_linear_joint'); % class name
% 
%     % The actuator with rigid gearbox and
%     % asymmetric coulomb and viscous friction 
%     jb.buildJoint('WMBig10k',...parameters
%         'rigid_gearbox',...  linear dynamics
%         {'coulomb_asym',...  nonlinear 
%          'viscous_asym'},... dynamics
%         [],...     static electric subsystem
%         'my_nonlinear_joint'); % class name
% 
%     % The same actuator with locked output 
%     %and a full linear electric dynamics
%     jb.buildJoint('WMBig10k', ... parameters
%         'output_fixed_rigid_gearbox', ... linear dynamics
%         {'coulomb_asym',...  nonlinear
%         'viscous_asym'}, ... dynamics
%         'electric_dyn',... electro-dynamics
%         'my_locked_joint');  % class name
% 
%     %% Instantiate joint models
%     % add build directory to search path
%     addpath(jb.buildDir) 
%     % declare joint objects
%     jointA = my_linear_joint;
%     jointB = my_nonlinear_joint;
%     jointC = my_locked_joint;
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also genericJoint.

% Copyright (C) 2016, by Joern Malzahn, Wesley Roozing
%
% This file is part of the Compliant Joint Toolbox (CJT).
%
% CJT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CJT is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
% License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CJT. If not, see <http://www.gnu.org/licenses/>.
%
% For more information on the toolbox and contact to the authors visit
% <https://github.com/geez0x1/CompliantJointToolbox>

classdef jointBuilder
    %JOINTBUILDER Builds joint model files
    
    properties
        buildDir = ['.',filesep,'build'];
        overwrite = 0;
        basePath;
        nonlinModelPath;
        linModelPath;
    end
    
    methods
        %__________________________________________________________________
        % Constructor
        function this = jointBuilder()
            
            % Gather information about the paths
            % 1) location of the joint builder
            pathstr = fileparts(which('jointBuilder.m'));
            this.basePath = cd(cd([pathstr, filesep, '..', filesep]));
            
            % 2) location of the build directory
            this.buildDir = [this.basePath, filesep, 'build'];
            
            % 3) location of the nonlinear models
            this.nonlinModelPath = [this.basePath, filesep, 'model', filesep, 'nonlinear'];
            
            % 4) location of the linear models
            this.linModelPath = [this.basePath, filesep, 'model', filesep, 'linear'];
        end
        
        %__________________________________________________________________
        % Remove all created files
        function purge(this)
            if (exist(this.buildDir, 'dir'))
                delete([this.buildDir filesep '*']);
            end
        end
        
        %__________________________________________________________________
        % Build a joint object file
        function className = buildJoint(this, paramName, modelName, nonlinearModelName, electricalDynamicsName, className)
        % buildJoint(paramName, modelName[, nonlinearModelName, electricalDynamicsName, className]) 
        % Builds a joint model file in the build/ directory based on
        % the supplied parameter and model names
        %
        % **paramName**: The name of an m-file defining a structure called 'param' with the fields defined in the
        % genericJoint class. The genericJoint class already specifies default values for all parameters. The 
        % parameter file therefore has to define those parameters deviating from the default. 
        %
        % **paramName** is a |REQUIRED| input parameter.
        %
        % **modelName**: This is the name of the linear model template to be used. The following linear model templates can
        % be selected here:
        %   * *full_dyn*: Three mass system with motor, gear and load inertia. Gearbox and sensor are compliant.
        %   * *full_dyn_no_friction*: As above, but without friction terms.
        %   * *output_fixed*: Two mass system with motor and gear inertia. The system output is locked.
        %   * *output_fixed_no_friction*: Same as above, but without friction terms.
        %   * *output_fixed_rigid_gearbox*: Single mass spring damper system with locked output.
        %   * *output_fixed_rigid_gearbox_no_friction*: As above, but without friction terms.
        %   * *rigid*: Single mass system without compliant elements.
        %   * *rigid_no_friction*: As above, but without friction terms.
        %   * *rigid_gearbox*: Two mass system with combined motor-gearbox ineartia.
        %   * *rigid_gearbox_no_friction*: As above, but without friction terms.
        %
        %
        % **modelName** is an |OPTIONAL| input parameter and defaults to *full_dyn*
        %
        % **nonlinearModelName**: Specifies the nonlinear dynamics added to the linear model. The parameter is of cell
        % type and can contain a list of multiple effects. Possible "ingredients" are:
        %
        %   * *coulomb*: Symmetric Coulomb friction, parameters d_cx are used.
        %   * *coulomb_asym*: Asymmetric Coulomb friction, parameters d_cx and d_cx_n are used.
        %   * *torque_ripple*: Harmonic torque ripple, parameters rip_a1, rip_a2 and rip_f are used.
        %   * *viscous_asym*: Asymmetric viscous friction, parameters d_x_n are used.
        %
        % **nonlinearModelName** is an |OPTIONAL| input parameter. Defaults to an empty list and purely linear dynamics
        % models.
        %
        % **electricalDynamicsName**: Specifies the dynamics to be considered for the electrical subsystem. Possible
        % choices are:
        %
        %  * *electric_dyn*: Linear electrical dynamics including armature resistance and inductance. May result in slow
        %  simulations!
        %  * *electric_dyn_zero_inductance*: Static model of the delectrical dynamics, inductance is ignored. 
        %
        % **electricalDynamicsName** is an |OPTIONAL| input parameter. Defaults to *electric_dyn_zero_inductance*.
        %
        % **className**: Allows to give the derived joint model a custom name.
        % **className**:  is an |OPTIONAL| input parameter. If no custom name is given, the method creates one from the
        % specified parameters and dynamics.
        

            
            % Check whether the params exist
            if (exist(paramName, 'file') ~= 2)
                error(['jointBuilder.buildJoint error: Parameters ''' paramName ''' do not exist!']);
            end

            % Evaluate parameters
            eval(paramName);
            
            % If model name is empty, use "full_dyn" by default
            if (~exist('modelName', 'var') || isempty(modelName))
                modelName = 'full_dyn'; 
            end
            % Otherwise, check whether the linear model exists
            % This relies on the model function name being equal to
            % the filename (which we require)
            if (~exist([this.linModelPath, filesep, modelName,'.m'], 'file'))
                error(['jointBuilder.buildJoint error: Linear model ''' modelName ''' do not exist!']);
            end
            
            % Process nonlinear model names to create the class name as
            % well as the nonlinear model property fo the class
            if (~exist('nonlinearModelName', 'var') || isempty(nonlinearModelName))
                % No nonlinear dynamics
                nonlinearModelName = '';
                nonlinearModelNameCellStr = '''''';
            else
                % One or multiple nonlinear dynamics components included
                
                % If the argument is a cell, we're dealing with multiple
                % nonlinear dynamics components
                if (iscell(nonlinearModelName) && length(nonlinearModelName) > 1)
                    % Check if all nonlinear models exist
                    for i=1:length(nonlinearModelName)
                        % Check whether it exists
                        % This relies on the model function name being equal to
                        % the filename (which we require)
                        if (~exist(nonlinearModelName{i}, 'file'))
                            error(['jointBuilder.buildJoint error: Nonlinear model ''' nonlinearModelName{i} ''' does not exist!']);
                        end
                    end
                    
                    % Multiple components
                    nonlinearModelNameStr = strjoin(nonlinearModelName, '_');
                    nonlinearModelNameCellStr = strcat('{''', strtrim(strjoin(nonlinearModelName, ''', ''')), '''}');
                else
                    % In case one nonlinear model was given, but as a cell
                    if (iscell(nonlinearModelName))
                        nonlinearModelName = nonlinearModelName{:};
                    end
                    
                    % Check whether it exists
                    % This relies on the model function name being equal to
                    % the filename (which we require)
                    if (~exist([this.nonlinModelPath, filesep, nonlinearModelName], 'file'))
                        error(['jointBuilder.buildJoint error: Nonlinear model ''' nonlinearModelName ''' does not exist!']);
                    end
                    
                    % Single component
                    nonlinearModelNameStr = nonlinearModelName;
                    nonlinearModelNameCellStr = ['''' nonlinearModelName ''''];
                    
                    % Turn nonlinearModelName into a cell so we can process
                    % them in the same way
                    nonlinearModelName = {nonlinearModelName};
                end
            end
            

            % Check whether the a linear model for electrical dynamics is
            % specified. If not use the zero inductance model.
            if (~exist('electricalDynamicsName', 'var')) || isempty(electricalDynamicsName)
                electricalDynamicsName = 'electric_dyn_zero_inductance';
            end
            % Check whether the linear model for electrical dynamics exists
            % This relies on the model function name being equal to
            % the filename (which we require)
            if (~exist(electricalDynamicsName, 'file'))
                error(['jointBuilder.buildJoint error: Electrical subsystem model ''' electricalDynamicsName ''' do not exist!']);
            end
            
            % If no used-defined name is given, create the class name
            if (~exist('className', 'var'))
                className = [paramName, '_', modelName];
                if ~isempty(nonlinearModelName)
                    className = [className, '_', nonlinearModelNameStr];
                end
                % if the electrical dynamics are not set to the default
                % zero inductance model, indicate that in the name:
                if ~strcmp(electricalDynamicsName,'electric_dyn_zero_inductance')
                    className = [className, '_', electricalDynamicsName];
                end
            end
            className(isspace(className)) = []; % Remove any spaces
            
            % Joint name is equal to the long class name since the name
            % character limit only applies to filenames
            jointName = className;
            
            % Fix long filenames
            classNameLong = className;
            if (length(className) > 63-2) % -2 for .m extension
                className = [className(1:63-2-3) '_xx']; % -3 for '_xx'
            end
            
            % Create build directory if necessary
            if ~exist(this.buildDir,'dir')
                mkdir(this.buildDir);
                addpath(this.buildDir);
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
            
            
            % Class definition
            fprintf(fid,'%% %s\n\n', classNameLong);
            fprintf(fid,'classdef %s < genericJoint\n\n', className);
            
            % Description
            % ...
            
            
            % Private properties
            fprintf(fid,'\tproperties (SetAccess = private)\n');
            fprintf(fid,'\tend\n\n');
            
            
            % Methods
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
            getNonlinDynStr = fileread('getNonlinearDynamics.m');
            getNonlinDynStr = regexprep(getNonlinDynStr, '^', '\t\t', 'emptymatch', 'lineanchors');
            if (~isempty(nonlinearModelName))
                % There are nonlinear terms
                fStr = '';
                for i=1:length(nonlinearModelName)
                   fStr = strcat(fStr, nonlinearModelName(i), '(obj, x)');
                   if (i < length(nonlinearModelName))
                       fStr = strcat(fStr, '+');
                   end
                end
                fStr = strcat(fStr, ';');
            else
                % No nonlinear terms
                fStr = 'zeros(size(x)); % No nonlinear dynamics!';
            end
            fprintf(fid, [getNonlinDynStr '\n\n'], char(fStr));
            
            % getElectricalDynamicsMatrices
            getDynStr = fileread('getElectricalDynamicsMatrices.m');
            getDynStr = regexprep(getDynStr, '^', '\t\t', 'emptymatch', 'lineanchors');
            fprintf(fid, [getDynStr '\n\n'], electricalDynamicsName);
            
            % End methods
            fprintf(fid,'\tend\n\n');
            
            % End class definition
            fprintf(fid,'end\n');
            
            % Close file
            fclose(fid);
            
            % Status
            disp(['Joint built and written to ' classFName]);
        end
        
    end
    
end

