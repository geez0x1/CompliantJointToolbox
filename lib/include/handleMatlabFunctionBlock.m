function [  ] = handleMatlabFunctionBlock( thisBlock, jointObj )
%HANDLEMATLABFUNCTIONBLOCK Handles the contents of the MatlabFunctionBlock
%inside the joint model block.
%
% [ ] = handleMatlabFunctionBlock( thisBlock, jointObj )
%
% Inputs::
% thisBlock Simulink block handle to the block containing the MatlabFunctionBlock
% jointObj Object defining the model to use.
%
%
% Notes::
%
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also .

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

% Get linear dynamics matrices
[A, B, C, I, D, K] = jointObj.getDynamicsMatrices();

% Other properties
Ts = jointObj.Ts;
k_t = jointObj.k_t;
k_b = jointObj.k_b;

% ---------------------------------
% Edit contents of nonLinDynamicsFcn block
% ---------------------------------
% thisBlock = gcb;
funName = 'nonLinDynamicsFcn';
fcnBlock = [thisBlock, '/', funName];
intBlock = [thisBlock, '/', 'Integrator'];
sumBlock = [thisBlock, '/', 'Sum'];
%refBlock = get_param(fcnBlock,'ReferenceBlock');

% Locate the property, which stores the source code text
% r = slroot;
%  scriptBlock = r.find('-isa','Stateflow.EMChart','path',fcnBlock);
r = slroot;
scriptBlock = r.find('-isa','Stateflow.EMChart','path',fcnBlock);
if isempty(scriptBlock)
    load_system('simulink');
    
    delete_line(thisBlock,['B_nonlinear/1'], ['Sum/1']);
    delete_line(thisBlock,[funName '/1'], ['B_nonlinear/1']);
    delete_line(thisBlock,[funName '/1'], ['nonLinearDynamicsScope/1']);
    delete_line(thisBlock,'Integrator/1',[funName '/1'])
    delete_block(fcnBlock)

    add_block('simulink/User-Defined Functions/MATLAB Function',fcnBlock,...
        'Position', [630 312 700 358],...
        'Orientation', 'left');
    add_line(thisBlock,[funName '/1'], ['B_nonlinear/1'],'AutoRouting','On');
    add_line(thisBlock,['B_nonlinear/1'], ['Sum/1'],'AutoRouting','On');
    add_line(thisBlock,[funName '/1'], ['nonLinearDynamicsScope/1'],'AutoRouting','On');
    add_line(thisBlock,'Integrator/1',[funName '/1'],'AutoRouting','On')
    

    scriptBlock = r.find('-isa','Stateflow.EMChart','path',fcnBlock);
    if isempty(scriptBlock)
        error(message('symbolic:sym:matlabFunctionBlock:CouldNotCreate', fcnBlock));
    end
end
if size(scriptBlock) > 1
    error(message('symbolic:sym:matlabFunctionBlock:AmbiguousBlock', fcnBlock));
end
if ~isa(scriptBlock,'Stateflow.EMChart')
    error(message('symbolic:sym:matlabFunctionBlock:InvalidBlock', fcnBlock));
end
 
 
% The strategy is to write a temporary text file to format the code.
% Then read the file again and paste the string into the source code
% property.
file = 'nonLinDynamicsFcnTMP.m';
fid = fopen(file,'wt+');
tmp = onCleanup(@()delete(file)); % file guard -> removes temporary file upon crash.

fprintf(fid,'%s\n',['function tau = nonLinDynamics(x)']);   % write signature
fprintf(fid,'%s\n\n',['%#codegen']);                            % indicate possibility of code generation

% Make parameters available inside the function
s = jointObj.getParams();   % Get params struct
fields = fieldnames(s)';    % Get fieldnames
N = numel(fields);          % Loop through all fields and print
for i=1:N
    val = s.(fields{i});
    if (isnumeric(val))
        % Numeric values
        if (length(val) > 1 || isempty(val))
            % Arrays (1-D)
            fprintf(fid, 'params.%s = [', fields{i});
            for j=1:length(val)
                fprintf(fid, '%f', val(j));
                if (j < length(val))
                    if (size(val,1) > size(val,2))
                        fprintf(fid, '; ');
                    else
                        fprintf(fid, ', ');
                    end
                end
            end
            fprintf(fid, '];\n');
        else
            % Simple floating point value
            fprintf(fid, 'params.%s = %f;\n', fields{i}, val);
        end
        
    elseif (iscell(val))
        % String cells
        fprintf(fid, '%%params.%s = {', fields{i});
        for j=1:length(val)
            fprintf(fid, '''%s''', s.(fields{i}){j});
            if (j < length(val))
                fprintf(fid, ', ');
            end
        end
        fprintf(fid, '}; %% Cell arrays are not supported for code generation.\n');
        
    else
        % Strings
        fprintf(fid, 'params.%s = ''%s'';\n', fields{i}, val);
        
    end
end

% Write the nonlinear model code
argStr = '(params,x)';
nonlinearModelName = jointObj.nonlinearModelName;
nonlinearModelName_multiple = 0;
if (iscell(nonlinearModelName) && length(nonlinearModelName) > 1)
    nonlinearModelName_multiple = 1;
end

fStr = '';
if (~isempty(nonlinearModelName) && nonlinearModelName_multiple)
    % Multiple nonlinear terms
    for i=1:length(nonlinearModelName)
        fStr = strcat(fStr, nonlinearModelName(i), argStr);
        if (i < length(nonlinearModelName))
            fStr = strcat(fStr, char(' + '));
        end
    end
elseif ~isempty(nonlinearModelName)
    % Single nonlinear term
    fStr = strcat(nonlinearModelName,argStr);
else
    % No nonlinear terms
    fStr = sprintf('%s',  'zeros(size(x))');
end
fStr = strcat('tau = ', fStr, ';');
fprintf(fid, '%s\n', char(fStr));

% End methods
fprintf(fid,'\tend\n\n');


contents = fileread(file);
scriptBlock.Script = sprintf('%s',contents);

fclose(fid);
clear('tmp');

outPortAddress = [fcnBlock, '/', 'tau'];
set_param(outPortAddress,'PortDimensions',num2str(size(A,1)));


end

