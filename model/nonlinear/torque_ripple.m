%COGGING Calculate torque ripple torques
%
% [tau] = torque_ripple(jointObj, x)
%
% jointObj is the instance of the joint class object for which this
% function has been called.
%
% Inputs::
%   x: state vector depending on the model type as
%     x = [q_m; q_g; q_l; q_m_dot; q_g_dot, q_l_dot'];  full_dyn
%     x = [q_g, q_l, q_g_dot, q_l_dot]'                 rigid_gearbox
%     x = [q_m, q_g, q_m_dot, q_g_dot]'                 output_fixed
%     x = [q_g, q_g_dot]'                               output_fixed_rigid_gearbox
%     x = [q_g, q_g_dot]'                               rigid
%
% Outputs::
%   tau: torque vector of appropriate size
%
% Notes::
%   Supported torque ripple types (indexes are values of rip_types):
%     1:  Position-dependent (e.g. magnetic cogging - 'detent' or 'no-current' torque)
%     2:  Position- and torque-dependent (e.g. Harmonic Drive torque ripple)
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

function [ tau ] = torque_ripple(jointObj, x)
    
    % Torque_ripple parameters
    types   = jointObj.rip_types;   % Torque ripple types
    a1      = jointObj.rip_a1;      % Cosine amplitudes [Nm]
    a2      = jointObj.rip_a2;      % Sine amplitudes [Nm]
    f       = jointObj.cog_f;       % Spatial frequencies [periods/revolution]
    omegas	= 2 * pi * f;           % Spatial frequencies [rad/revolution]
    
    % Preallocate coefficient vector
    c = zeros(size(x));
    
    % Calculate torque for all ripple sources
    n_types     = length(types);
    taus        = zeros([n_types 1]);
    for i=1:n_types
        if (types(i) == 1)          % Position-dependent
            taus(i) = a1(i) * cos(omegas(i) * x(1)) + a2(i) * sin(omegas(i) * x(1));
        elseif (types(i) == 2)      % Position- and torque-dependent
            taus(i) = 0;
        else
            error(['ERR: Invalid torque ripple type specified: types(' num2str(i) ') = ' types(i)]);
        end
    end
    
    % Build coefficient vector
    if (strcmp(jointObj.modelName, 'full_dyn'))
        c = [   0;
                0;
                0;
                sum(taus);
                0;
                0       	];
        
    elseif (strcmp(jointObj.modelName, 'rigid_gearbox'))
        c = [   0;
                0;
                sum(taus);
                0           ];
        
    elseif (strcmp(jointObj.modelName, 'output_fixed'))
        c = [   0;
                0;
                sum(taus);
                0       	];
        
    elseif (strcmp(jointObj.modelName, 'output_fixed_rigid_gearbox'))
        c = [   0;
                sum(taus)   ];
        
    elseif (strcmp(jointObj.modelName, 'rigid'))
        c = [   0;
                sum(taus)   ];
        
    end

    % Calculate cogging torques
    tau = c;

end

