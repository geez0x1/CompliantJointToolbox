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
        function className = buildJoint(this, paramName, modelName, nonlinearModelName, className)
            % buildJoint(paramName, modelName, nonlinearModelName)
            % Builds a joint model file in the build/ directory based on
            % the supplied parameter and model names

            % Eval parameters
            eval(paramName);
            
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
            if (~exist('className', 'var'))
                className = [paramName '_' modelName];
                if ~isempty(nonlinearModelName)
                    className = [className '_' nonlinearModelNameStr];
                end
            end
            className(isspace(className)) = []; % Remove any spaces
            
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

