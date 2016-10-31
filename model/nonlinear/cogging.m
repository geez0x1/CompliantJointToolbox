%COGGING Calculate magnetic cogging torque
%
% [tau] = cogging(jointObj, x)
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
%   tau: cogging torque vector of appropriate size
%
% Notes::
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also 

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

function [ tau ] = cogging(jointObj, x)
    
    % Cogging parameters
    a1      = jointObj.cog_a1;      % Cosine amplitude [Nm]
    a2      = jointObj.cog_a2;      % Sine amplitude [Nm]
    f       = jointObj.cog_f;       % Spatial frequency [periods/revolution]
    omega	= 2 * pi * f;           % Spatial frequency [rad/revolution]
    
    % Preallocate coefficient vector
    c = zeros(size(x));
    
    % Build coefficient vector
    if (strcmp(jointObj.modelName, 'full_dyn'))
        c = [   0;
                0;
                0;
                a1 * cos(omega * x(1)) + a2 * sin(omega * x(1));
                0;
                0                                           ];
        
    elseif (strcmp(jointObj.modelName, 'rigid_gearbox'))
        c = [   0;
                0;
                a1 * cos(omega * x(1)) + a2 * sin(omega * x(1));
                0                                           ];
        
    elseif (strcmp(jointObj.modelName, 'output_fixed'))
        c = [   0;
                0;
                a1 * cos(omega * x(1)) + a2 * sin(omega * x(1));
                0                                           ];
        
    elseif (strcmp(jointObj.modelName, 'output_fixed_rigid_gearbox'))
        c = [   0;
                a1 * cos(omega * x(1)) + a2 * sin(omega * x(1))     ];
        
    elseif (strcmp(jointObj.modelName, 'rigid'))
        c = [   0;
                a1 * cos(omega * x(1)) + a2 * sin(omega * x(1)) 	];
        
    end

    % Calculate cogging torques
    tau = c;

end

