% GET_CONTROLLED_CLOSED_LOOP Outputs closed loop transfer functions.
%
%   [ P, G, H, Kd_opt ] = get_controlled_closed_loop( jointObj, Kp, Ki, Kd, N [, pid_form, outputIdx, ff_comp_switch] )
%
%   Gets a plant P, PD controller G, PD+(FF or compensation) closed loop
%   H transfer functions, and the PD derivative gain Kd_opt that critically
%   damps the closed-loop poles for fixed-output force control.
%
% Inputs::
%   jointObj: Joint object
%   Kp, Ki, Kd, N: PID gains and derivative cut-off frequency (if Kd=-1, the
%                  resulting controller D-gain will be set to Kd_opt).
%   pid_form: Flag that determines whether PID controller is constructed in
%             product (ideal) or summation (parallel) form
%   outputIdx: Controlled/observed plant output (default: 7, torque)
%   ff_comp_switch: Feed-forward/compensation (1=Compensation, 2=Feed-forward)
%
% Outputs::
%   P: Plant transfer function
%   G: PID Controller transfer function
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

function [ P, G, H, Kd_opt ] = get_controlled_closed_loop(jointObj, Kp, Ki, Kd, N, pid_form, outputIdx, ff_comp_switch)
    %% Default arguments
    if (~exist('pid_form', 'var') || isequal(pid_form,[]))
        pid_form = 'ideal';     % Ideal PID form (series) by default
    end
    if (~exist('outputIdx', 'var') || isequal(outputIdx,[]))
        outputIdx = 7;          % Controlled/observed plant output (default: 7, torque)
    end
    if (~exist('ff_comp_switch', 'var') || isequal(ff_comp_switch,[]))
        ff_comp_switch = 1;     % Feed-forward/compensation
                                % (1=Compensation (default), 2=Feed-forward)
    end
    
    
    %% Get variables
    
    % Control/system parameters
    k_b     = jointObj.k_b;     % Torsion bar stiffness [Nm/rad]
    n       = jointObj.n;       % Gearbox transmission ratio []
    k_t     = jointObj.k_t;     % Torque constant [Nm/A]
    I_m     = jointObj.I_m;     % Motor rotor inertia [kg m^2]   
    I_g     = jointObj.I_g;     % Gear inertia [kg m^2]
    d_m     = jointObj.d_m;     % Motor Damping [Nms/rad]
    d_g     = jointObj.d_g;     % Gearbox damping [Nms/rad]
    d_gl    = jointObj.d_gl;    % Torsion bar internal damping [Nms/rad]
    
    % Optimal (specified damping ratio of poles) gains
    zeta = 1.0; % Critical damping
    Kd_opt =    (   2 * (I_m + I_g) * ...
                    sqrt( ...
                        (k_b * Kp + k_b) / ...
                        (I_m + I_g) ...
                    ) * zeta ...
                    - (d_m + d_g + d_gl) ...
                ) / k_b;
    
    % Set derivative gain Kd to Kd_opt if set to -1
    if (Kd == -1)
        Kd = Kd_opt;
    end
    
    
    %% Get state-space system with current input and specified output
    sys         = jointObj.getStateSpace();
    sys         = ss(sys.A, sys.B(:,1), sys.C(outputIdx,:), sys.D(outputIdx,1));
    
    
    %% Build closed-loop system
    
    % PID controller
    if (strcmpi(pid_form, 'ideal'))
        G = pidstd(Kp, 1/(Kp*Ki), Kd/Kp, N);    % ideal PID controller
    elseif (strcmpi(pid_form, 'parallel'))
        G = pid(Kp, Ki, Kd, 1/N);               % parallel PID controller
    else
        error('Invalid PID form');
    end
    G.u    = 'e';
    G.y    = 'pid_u';

    
    % See whether we need to build a closed-loop system with compensation or
    % feed-forward for the spring force
    if (ff_comp_switch == 1)            % Compensation
        % Compensation term
        C2      = tf(1);
        
        % Closed-loop transfer function
        P       = tf(sys);
        P2      = feedback(P, C2, +1);  % Positive feedback for the comp/FF
        Pf      = feedback(P2 * G, 1);
        
    elseif (ff_comp_switch == 2)        % Feed-forward
        % Feed-forward term
        C2      = tf(1);
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
        P       = tf(sys);
        Pf      = feedback(P * G, 1);
        
    else
        error('Error: ff_comp_switch not set to 1, 2 or 3.');
    end
    
    % The closed loop H is equal to Pf
    H = tf(Pf);
    
end

