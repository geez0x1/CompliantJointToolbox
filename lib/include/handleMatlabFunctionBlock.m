function [  ] = handleMatlabFunctionBlock( thisBlock, jointObj )
%HANDLEMATLABFUNCTIONBLOCK Summary of this function goes here
%   Detailed explanation goes here


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
    
    delete_line(thisBlock,[funName '/1'], ['Sum/1']);
    delete_line(thisBlock,[funName '/1'], ['nonLinearDynamicsScope/1']);
    delete_line(thisBlock,'Integrator/1',[funName '/1'])
    delete_block(fcnBlock)

    add_block('simulink/User-Defined Functions/MATLAB Function',fcnBlock,...
        'Position', [630 312 700 358],...
        'Orientation', 'left');
    add_line(thisBlock,[funName '/1'], ['Sum/1'],'AutoRouting','On');
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
        fprintf(fid, 'params.%s = %f;\n', fields{i}, s.(fields{i}));
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
        fprintf(fid, 'params.%s = ''%s'';\n', fields{i}, s.(fields{i}));
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

