% GET_CONTROLLED_CLOSED_LOOP Outputs closed loop transfer functions.
%
%   [ P, G, H, Kd_opt ] = get_controlled_closed_loop( jointName, Kp, Ki, Kd, N [, ff_comp_switch] )
%
%   Gets a plant P, PD controller G, PD+(FF or compensation) closed loop
%   H transfer functions, and the PD derivative gain Kd_opt that critically
%   damps the closed-loop poles for fixed-output force control.
%
% Inputs::
%   jointName: Joint object class name Kp, Ki, Kd, N: PID gains and 
%              derivative cut-off frequency (if Kd=-1, the
%              resulting PD D-gain will be set to Kd_opt).
%   ff_comp_switch: Feed-forward/compensation (1=Compensation, 2=Feed-forward)
%
% Outputs::
%   P: Plant transfer function
%   G: Controller transfer function
%   H: Closed-loop transfer function
%   Kd_opt: Optimal PD derivative gain that critically damps the closed-loop poles for
%           fixed-output force control.
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
% See also getObserver, getKalman.

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


function [ P, G, H, Kd_opt ] = get_controlled_closed_loop( jointName, Kp, Ki, Kd, N, ff_comp_switch )
    %% Default arguments
    
    if (~exist('ff_comp_switch', 'var'))
        ff_comp_switch = 1;     % Feed-forward/compensation
                                % (1=Compensation (default), 2=Feed-forward)
    end

    
    %% Get variables

    % Get joint object
    j = eval(jointName);

    
    %% Control/system parameters
    k_b    	= j.k_b;                % Torsion bar stiffness [Nm/rad]
    n       = j.n;                  % Gearbox transmission ratio []
    k_t     = j.k_t;                % Torque constant [Nm/A]
    
    % Optimal (specified damping ratio of poles) gains
    zeta = 1.0; % Critical damping
    Kd_opt = 	(   2*(j.I_m + j.I_g) * ...
                    sqrt( ...
                        (j.k_b * j.k_t * Kp * j.n + j.k_b) / ...
                        (j.I_m + j.I_g) ...
                    ) * zeta ...
                    - (j.d_m + j.d_g + j.d_gb) ...
                ) / (j.k_b * j.k_t * Kp * j.n);
    
    % Set derivative gain to Kd_opt if set to -1
    if (Kd == -1)
        Kd = Kd_opt;
    end
    
  
    %% Get state-space system with torque output
    sys         = j.getStateSpace();
    sys         = ss(sys.A, sys.B, sys.C(2,:), 0);
    sys         = k_b * sys; % Multiply by k_b to get torque output
    
    
    %% Build closed-loop system

    % Controllers
    G    	= pid(Kp, Kp*Ki, Kp*Kd, 1/N); % PID controller
    G.u  	= 'e';
    G.y 	= 'pid_u';

    % See whether we need to build a closed-loop system with compensation or
    % feed-forward for the spring force
    if (ff_comp_switch == 1)        % Compensation
        % Compensation term
        C2 = 1 / (n * k_t);

        % Closed-loop transfer function
        P	= tf(sys);
        P2	= feedback(P, C2, +1); % Positive feedback
        Pf	= feedback(P2 * G, 1);

    elseif (ff_comp_switch == 2)    % Feed-forward
        % Feed-forward term
        C2      = tf(1 / (n * k_t));
        C2.u    = 'r';
        C2.y    = 'ff_u';

        % Closed-loop transfer function
        P       = tf(sys);
        P.u     = 'u';
        P.y     = 'y';
        sum_e   = sumblk('e = r - y');
        sum_ff  = sumblk('u = pid_u + ff_u');
        Pf      = connect(P, G, C2, sum_e, sum_ff, 'r', 'y');

    elseif (ff_comp_switch == 3)    % Spring force compensation off
        % Closed-loop transfer function
        P	= tf(sys);
        Pf	= feedback(P * G, 1);

    else
        error('Error: ff_comp_switch not set to 1, 2 or 3.');
        return;
    end


    %% Pc = Pf
    H = tf(Pf);

    
end

