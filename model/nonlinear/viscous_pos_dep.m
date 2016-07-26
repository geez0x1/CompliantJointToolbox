%VISCOUS_POS_DEP Calculate position-dependent viscous friction
%
% [tau] = jointObj.viscous_pos_dep(x)
%
% jointObj is the instance of the joint class object for which this
% function has been called.
%
% Inputs::
%   x: state vector depending on the model type as
%     x = [q_m; q_g; q_b; q_m_dot; q_g_dot, q_b_dot'];  full_dyn
%     x = [q_g, q_b, q_g_dot, q_b_dot]'                 rigid_gearbox
%     x = [q_m, q_g, q_m_dot, q_g_dot]'                 output_fixed
%     x = [q_g, q_g_dot]'                               output_fixed_rigid_gearbox
%     x = [q_g, q_g_dot]'                               rigid
%
% Outputs::
%   tau: friction torque
%
% Notes::
%  THIS FILE IS A PROOF OF CONCEPT AND NEEDS TO BE FINISHED
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also full_dyn, coulomb, viscous_asym.

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

function [ tau ] = viscous_pos_dep(obj, x)
    
    % Temp
    a = 1e-6;   % Amplitude
    d = 1/2;    % Frequency
    p = pi;     % Phase
    
    c = zeros(size(x));
    
    % Build coefficient vector
    if (strcmp(obj.modelName, 'full_dyn'))
        c = [0, 0, 0, ...
                a * sin(x(1)/d + p), ...
                a * sin(x(2)/d + p), ...
                a * sin(x(3)/d + p)       ]';
        
    elseif (strcmp(obj.modelName, 'rigid_gearbox'))
        c = [0, 0, ...
                a * sin(x(1)*d + p) + a * sin(obj.n*x(1)*d + p) * obj.n, ...
                a * sin(x(2)*d + p)       ]';
        
    elseif (strcmp(obj.modelName, 'output_fixed'))
        c = [0, 0, ...
                a * sin(x(1)/d + p), ...
                a * sin(x(2)/d + p)       ]';
        
    elseif (strcmp(obj.modelName, 'output_fixed_rigid_gearbox'))
        c = [0, a * sin(x(1)*d + p) + a * sin(obj.n*x(1)*d + p) * obj.n]';
        
    elseif (strcmp(obj.modelName, 'rigid'))
        c = a * sin(x(1)*d + p) + a * sin(obj.n*x(1)*d + p) * obj.n;
        
    end

    % Calculate position-dependent viscous friction torques
    tau = -c .* x;

end

