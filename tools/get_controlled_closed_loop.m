% GET_CONTROLLED_CLOSED_LOOP Outputs closed loop transfer functions.
%
%   [ P, H, Kd_opt ] = get_controlled_closed_loop( jointObj, Kp, Ki, Kd, N [, pid_form, outputIdx, ff_comp_switch, derivative_select, zeta] )
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
%   ff_comp_switch: Feed-forward/compensation ('comp' (compensation, default), 'ff' (feed-forward), or 'off')
%   derivative_select: Derivative action on error or output ('error' (default) or 'output')
%   zeta: Pole damping to use for Kd_opt (default: 1.0)
%
% Outputs::
%   P: Plant transfer function
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

function [ P, H, Kd_opt ] = get_controlled_closed_loop(jointObj, Kp, Ki, Kd, N, pid_form, outputIdx, ff_comp_switch, derivative_select, zeta)
    %% Default arguments
    if (~exist('pid_form', 'var') || isequal(pid_form,[]))
        pid_form = 'ideal';             % Ideal PID form (series) by default
    end
    if (~exist('outputIdx', 'var') || isequal(outputIdx,[]))
        outputIdx = 7;                  % Controlled/observed plant output (default: 7, torque)
    end
    if (~exist('ff_comp_switch', 'var') || isequal(ff_comp_switch,[]))
        ff_comp_switch = 'comp';        % Feed-forward/compensation
                                        % ('comp' (compensation, default), 'ff' (feed-forward), or 'off')
    end
    if (~exist('derivative_select', 'var') || isequal(derivative_select,[]))
        derivative_select = 'error';    % Derivative action on error or output
                                        % ('error' (default) or 'output')
    end
    if (~exist('zeta', 'var') || isequal(zeta,[]))
        zeta = 1.0;                     % Pole damping to use for Kd_opt (default: 1.0)
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
    Kd_opt =    (   2 * (I_m + I_g) * ...
                    sqrt( ...
                        (k_b * Kp + k_b) / ...
                        (I_m + I_g) ...
                    ) * zeta ...
                    - (d_m + d_g + d_gl) ...
                ) / k_b;
    
    % Set derivative gain Kd to Kd_opt if set to -1
    if (Kd == -1)
        % The above optimal D-gain was computed for parallel PID form.
        % Convert it for ideal PID form.
        if (strcmpi(pid_form, 'ideal'))
            Kd = Kd_opt/Kp;
        else
            Kd = Kd_opt;
        end
    end
    
    
    %% Get state-space system with current input and specified output
    sys         = jointObj.getStateSpace();
    sys         = ss(sys.A, sys.B(:,1), sys.C(outputIdx,:), sys.D(outputIdx,1));
    
    
    %% Build closed-loop system
    
    % PI controller
    if (strcmpi(pid_form, 'ideal'))
        G1 = pidstd(Kp, 1/Ki, 0, N);   % ideal PI controller
    elseif (strcmpi(pid_form, 'parallel'))
        G1 = pid(Kp, Ki, 0, 1/N);      % parallel PI controller
    else
        error('Invalid PID form');
    end
    G1.u    = 'e';
    G1.y    = 'pid_pi_u';
    
    % Derivative action
    % Check the PID form to set the D-gain correctly
    if (strcmpi(pid_form, 'ideal'))
        Kd_D = Kp * Kd; % Multiply by Kp for 'ideal' form
    else
        Kd_D = Kd; % Parallel form is independent
    end
    if (strcmpi(derivative_select, 'error'))
        G2      = pid(0, 0, Kd_D, 1/N); % D controller
        G2.u    = 'e'; % derivative action on error
    elseif (strcmpi(derivative_select, 'output'))
        G2      = pid(0, 0, -Kd_D, 1/N); % D controller (negative gain for sign of y)
        G2.u    = 'y'; % derivative action on output
    else
        error('Invalid value for derivative_select.');
    end
    G2.y    = 'pid_d_u';

    
    % See whether we need to build a closed-loop system with compensation or
    % feed-forward for the spring force
    if (strcmpi(ff_comp_switch, 'comp'))    % Compensation
        % Compensation term
        C2      = tf(1);
        C2.u    = 'y';
        C2.y    = 'comp_u';
        
        % Closed-loop transfer function
        P       = tf(sys);
        P.u     = 'u';
        P.y     = 'y';
        sum_e   = sumblk('e = r - y');
        sum_pid = sumblk('pid_u = pid_pi_u + pid_d_u');
        sum_comp= sumblk('u = pid_u + comp_u');
        H       = connect(P, G1, G2, C2, sum_e, sum_pid, sum_comp, 'r', 'y');
        
    elseif (strcmpi(ff_comp_switch, 'ff'))  % Feed-forward
        % Feed-forward term
        C2      = tf(1);
        C2.u    = 'r';
        C2.y    = 'ff_u';
        
        % Closed-loop transfer function
        P       = tf(sys);
        P.u     = 'u';
        P.y     = 'y';
        sum_e   = sumblk('e = r - y');
        sum_pid = sumblk('pid_u = pid_pi_u + pid_d_u');
        sum_ff  = sumblk('u = pid_u + ff_u');
        H       = connect(P, G1, G2, C2, sum_e, sum_pid, sum_ff, 'r', 'y');
        
    elseif (strcmpi(ff_comp_switch, 'off')) % Spring force compensation off
        % Closed-loop transfer function
        P       = tf(sys);
        P.u     = 'u';
        P.y     = 'y';
        sum_e   = sumblk('e = r - y');
        sum_pid = sumblk('u = pid_pi_u + pid_d_u');
        H       = connect(P, G1, G2, sum_e, sum_pid, 'r', 'y');
        
    else
        error('ff_comp_switch not set to ''comp'', ''ff'', or ''off''.');
    end
    
    % Cast to transfer function
    H = tf(H);
    
end

