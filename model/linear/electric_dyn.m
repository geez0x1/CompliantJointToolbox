%ELECTRIC_DYN Get linear dynamics matrices for the electrical subsystem
%
% [A, B, C, D] = jointObj.electric_dyn
%
% jointObj is the instance of the joint class object for which this
% function has been called.
%
%
% Outputs::
%   A:   System matrix
%   B:   Input matrix
%   C:   Output matrix
%   D:   Direct Feedthrough matrix
%
% Notes::
% When this model is used for simulation, it may require very small
% sampling times that prolong the computation time for the simulation.
% In the case where motor current transients are not important, they may be
% neglected by using the electric_dyn_zero_inductance model instead.
% 
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also electric_dyn_zero_inductance, full_dyn, output_fixed.

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


function [A, B, C, D] = electric_dyn(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [i]', where i is the motor current.

    R = obj.r;
    L = obj.x;
    k_t = obj.k_t * obj.n; % The geared torque constant.

    % State-space matrix
    A = -R/L;
    
    % u = [v, dot{q_m}] input voltage v and motor velocity dot{q_m} causing
    % back emf
    B = [1/L, -k_t];
    
    % Output
    C = 1;
    
    % Direct Feedthrough
    D = [0, 0];

end
