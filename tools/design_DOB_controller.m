% DESIGN_DOB_CONTROLLER Designs all parameters of a DOB based torque control
% scheme.
%
% [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = design_DOB_controller(jointObj, Kp, Ki, Kd, N [, pid_form, outputIdx, ff_comp_switch, f_c_FF, f_c_DOB])
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
% filter PQ_td, and feed-forward plant in version + filter PQ_ff.
%
% Inputs::
%   jointObj: Joint object
%   Kp: Inner control loop proportional gain
%   Ki: Inner control loop integral gain
%   Kd: Inner control loop derivative gain
%   N: PID derivative filter
%   pid_form: Flag that determines whether PID controller is constructed in
%             product or summation form
%   outputIdx: Controlled/observed plant output (default: 7, torque)
%   ff_comp_switch: Flag that determins whether torque feedforward is
%                   active or not (default: true)
%   f_c_FF: Feed-forward cutoff frequency [Hz] (default: 40)
%   f_c_DOB: DOB cutoff frequency [Hz] (default: 60)
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


function [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = design_DOB_controller(jointObj, Kp, Ki, Kd, N, pid_form, outputIdx, ff_comp_switch, f_c_FF, f_c_DOB)
    %% Default parameters
    if (~exist('pid_form', 'var'))
        pid_form = 'ideal';     % Ideal PID form (series) by default
    end
    if (~exist('outputIdx', 'var'))
        outputIdx = 7;          % Controlled/observed plant output (default: 7, torque)
    end
    if (~exist('ff_comp_switch', 'var'))
        ff_comp_switch = 1;     % Feed-forward/compensation
                                % (1=Compensation (default), 2=Feed-forward)
    end
    if (~exist('f_c_FF', 'var'))
        f_c_FF	= 40;           % Feed-forward cutoff frequency [Hz]
    end
    if (~exist('f_c_DOB', 'var'))
        f_c_DOB	= 60;           % DOB cutoff frequency [Hz]
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
    [~, ~, Pc, ~] =  get_controlled_closed_loop(jointObj, Kp, Ki, Kd, N, pid_form, outputIdx, ff_comp_switch);


    %% Design low-pass Butterworth filters

    % Q_td
    [a, b] = butter(order(Pc), omega_c_DOB, 's');
    Q_td = tf(a,b);

    % Q_ff
    [a, b] = butter(order(Pc), omega_c_FF, 's');
    Q_ff = tf(a,b);

    % Pc^-1 * Q_td
    PQ_td = inv(Pc) * Q_td;

    % Pc^-1 * Q_ff
    PQ_ff = inv(Pc) * Q_ff;
    
    
    %% Show Bode plots of results
%     figure(2); clf; hold on;
%     bode(Gf,        bodeOpt);
%     bode(Pc,        bodeOpt);
%     bode(inv(Pc),	bodeOpt);
%     bode(Q_td,      bodeOpt);
%     bode(Q_ff,      bodeOpt);
%     bode(PQ_td,     bodeOpt);
%     bode(PQ_ff,     bodeOpt);
%     xlim([0.1 100]);
%     grid on;
%     legend('Gf', 'P_c', 'P_c^{-1}', 'Q_{td}', 'Q_{ff}', 'PQ_{td}', 'PQ_{ff}');


    %% Optionally save results
    
    % Disabled as this causes problems when calling automatically from
    % Simulink masks
%     fname = 'design_DOB_controller_results.mat';
%     if confirm(['Do you want to save the results to ' fname ' [y/N]?'], 0)
%         save(fname, 'Pc', 'Q_td', 'Q_ff', 'PQ_td', 'PQ_ff');
%         disp(['Data saved to ' fname]);
%     end
    
end
