% DESIGN_DOB_CONTROLLER Designs all parameters of a DOB based torque control
% scheme, with the DOB as an OUTER loop (closed-loop DOB).
%
% [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = design_DOB_controller(jointObj, Kp, Ki, Kd, N [, pid_form, outputIdx, ff_comp_switch, derivative_select, f_c_FF, f_c_DOB, DOB_order])
%
% This function calculates the approximated closed-loop transfer function
% Pc, low-pass Q-filters, and the inverted models for a DOB with premulti-
% plication control scheme. This is the standalone code, that does not
% require measurement data from a Simulink model to design the model.
% Instead, it uses a linear state-space model to estimate the transfer
% function of the plant, thus the results are fundamentally different from
% the experimental version.
% This function returns the approximated closed-loop transfer function
% Pc, DOB filter Q_td, feed-forward filter Q_ff, DOB plant inversion +
% filter PQ_td, and feed-forward plant inversion + filter PQ_ff. If plant
% inversion is not intended to be used PQ_ff can simply be discarded.
%
% Inputs::
%   jointObj: Joint object
%   Kp: Inner control loop proportional gain
%   Ki: Inner control loop integral gain
%   Kd: Inner control loop derivative gain - set to -1 for critically
%       damped poles for the 2nd-order system (see
%       get_controlled_closed_loop)
%   N: PID derivative filter
%   pid_form: Flag that determines whether PID controller is constructed in
%             product or summation form (default: ideal/series form)
%   outputIdx: Controlled/observed plant output (default: 7, torque)
%   ff_comp_switch: Feed-forward/compensation (1=Compensation, 2=Feed-forward)
%   derivative_select: Derivative action on error or output (1 = error, 2 = output)
%   f_c_FF: Feed-forward cutoff frequency [Hz] (default: 40)
%   f_c_DOB: DOB cutoff frequency [Hz] (default: 60)
%   DOB_order: DOB order (>= relative order of plant+controller)
%
% Outputs::
%   Pc: Estimated closed-loop transfer function
%   Q_td: DOB Q-filter
%   Q_ff: Feed-forward Q-filter
%   PQ_td: Inverted plant + DOB Q-filter
%   PQ_ff: Inverted plant + Feed-forward Q-filter
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
% See also getLinearDOB, getKalman.

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


function [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = design_DOB_controller(jointObj, Kp, Ki, Kd, N, pid_form, outputIdx, ff_comp_switch, derivative_select, f_c_FF, f_c_DOB, DOB_order)
    %% Default parameters
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
    if (~exist('derivative_select', 'var') || isequal(derivative_select,[]))
        derivative_select = 1;  % Derivative action on error or output
                                % (1 = error (default), 2 = output)
    end
    if (~exist('f_c_FF', 'var') || isequal(f_c_FF,[]))
        f_c_FF  = 40;           % Feed-forward cutoff frequency [Hz]
    end
    if (~exist('f_c_DOB', 'var') || isequal(f_c_DOB,[]))
        f_c_DOB = 60;           % DOB cutoff frequency [Hz]
    end
    if (~exist('DOB_order', 'var') || isequal(DOB_order,[]))
        DOB_order = 0;          % DOB order []
    end

    % Bode options
    bodeOpt = bodeoptions;
    bodeOpt.FreqUnits = 'Hz';
    

    %% Get variables

    % Cut-off frequencies
    omega_c_FF  = 2 * pi * f_c_FF;  % Feed-forward (model inv) LPF cutoff frequency [rad/s]
    omega_c_DOB = 2 * pi * f_c_DOB; % DOB LPF cutoff frequency [rad/s]
    

    %% Build closed-loop system

    % Get controlled closed loop dynamics
    [~, ~, Pc, ~] = get_controlled_closed_loop(jointObj, Kp, Ki, Kd, N, pid_form, outputIdx, ff_comp_switch, derivative_select);

    % Check selected DOB order
    % If it was set to zero before (no argument given), set it to the
    % relative order of P.
    if (DOB_order == 0)
        DOB_order = relativeOrder(Pc);
    end
    if (DOB_order < relativeOrder(Pc))
        error('DOB order needs to be equal to or larger than the relative order of the plant+controller Pc.');
    end

    
    %% Design low-pass Butterworth filters

    % Q_td
    [a, b] = butter(DOB_order, omega_c_DOB, 's');
    Q_td = tf(a,b);

    % Q_ff
    [a, b] = butter(DOB_order, omega_c_FF, 's');
    Q_ff = tf(a,b);
    
    % Design DOB TFs by inverting the controlled plant dynamics

    % Pc^-1 * Q_td
    PQ_td = inv(Pc) * Q_td;

    % Pc^-1 * Q_ff
    PQ_ff = inv(Pc) * Q_ff;
    
end
