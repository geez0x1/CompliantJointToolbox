%ELECTRIC_DYN_DQ Get compute nonlinear dynamics for the electrical
%subsystem with a d-q axis current model
%
% [x_dot, tau_m] = electric_dyn_dq(jointObj, x, u)
%
% Inputs::
%   jointObj: instance of the joint class object for which this function
%             has been called.
%   x:        Current state vector x = [i_q, i_d]'
%   u:        Input voltages v_q, v_d, and motor velocity q_m_dot, given as
%             u = [v_q, v_d, q_m_dot]
%
% Outputs::
%   x_dot:    State derivative
%   tau_m:    Generated motor torque
%
%
% Notes::
% When this model is used for simulation, it may require very small
% sampling times that prolong the computation time for the simulation. If
% the d-axis is not important one may refer to electric_dyn instead, and if
% the electrical transients may be neglected, refer to the
% electric_dyn_zero_inductance model.
% 
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also electric_dyn, electric_dyn_zero_inductance.

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


function [x_dot, tau_m] = electric_dyn_dq(jointObj, x, u)
    
    % Hack for initialisation - problems with dimensioning
    if (length(x) ~= 2 || length(u) ~= 3)
        x_dot = [0;0];
        tau_m = 0;
        return;
    end

    % The computations below assume a state vector definition according to:
    % x = [i_q, i_d]', where i is the motor current.
    i_q         = x(1);
    i_d         = x(2);
    
    % u = [v_q, v_d, q_m_dot] input voltages v_q, v_d, and motor velocity
    % q_m_dot causing back-EMF
    v_q         = u(1);
    v_d         = u(2);
    q_m_dot     = u(3);

    % Get parameters
    r           = jointObj.r;   % Winding resistance [Ohm]
    L_q         = jointObj.x;   % Winding inductance [H] % Equal q-d inductances for SPMSMs
    L_d         = jointObj.x;   % Winding inductance [H] % Equal q-d inductances for SPMSMs
    k_t         = jointObj.k_t; % Torque constant [Nm/A]
    n           = jointObj.n;   % Gearbox transmission ratio []
    p           = jointObj.p;   % Number of pole pairs
    v_0         = jointObj.v_0; % Supply voltage [V]
    
    % Calculate dependent parameters
    % Here p denotes the number of pole pairs
    lambda_PM   = (2/3) * k_t / p;
    
    % Calculate dependent variables
    omega_m     = n * q_m_dot; % rotor angular velocity, reflected across gearbox
    
    % Calculate state derivative
    % x = [i_q; i_d];
    x_dot = [   (1/L_q) * v_q - (r/L_q) * i_q - (L_d/L_q) * p * omega_m * i_d - (lambda_PM * p * omega_m)/L_q;    ...
                (1/L_d) * v_d - (r/L_d) * i_d + (L_q/L_d) * p * omega_m * i_q;                                    ];
            
    % Calculate output
    tau_m = (3/2) * p * (lambda_PM * i_q + (L_d - L_q) * i_d * i_q);
    tau_m = n * tau_m; % Reflect across gearbox as is expected by the jointModel

end
