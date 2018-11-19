%NONONLINEARYDNAMICS Function that returns empty tau and y of appropriate
%size, used for when a model has no nonlinear dynamics
%
% [ tau, y ] = no_nonlinear_dynamics(jointObj, x)
%
% jointObj is the instance of the joint class object for which this
% function has been called.
%
%
% Inputs::
%   x: state vector depending on the model type as
%     x = [q_m; q_g; q_l; q_m_dot; q_g_dot, q_l_dot'];  full_dyn
%     x = [q_g, q_l, q_g_dot, q_l_dot]'                 rigid_gearbox
%     x = [q_m, q_g, q_l, q_m_dot, q_g_dot]'            output_fixed
%     x = [q_g, q_l, q_g_dot]'                          output_fixed_rigid_gearbox
%     x = [q_g, q_g_dot]'                               rigid
%
% Outputs::
%   tau: friction torque (all zero)
%   y:   output components (all zero)
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

function [ tau, y ] = no_nonlinear_dynamics(obj, x)
    
    % No dynamic torques
    tau = zeros(size(x));
    
    % No nonlinear output components
    y = zeros(7,1);

end

