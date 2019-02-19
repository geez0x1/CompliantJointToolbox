% MASK_PID_FFCOMP_OUTERLOOPDOB Mask code for the PID+FF/comp with outer-loop DOB block
%
% Calculate the approximated closed-loop transfer function Pc and
% filters Q_* etc for an outer-loop DOB, that is, a DOB that is closed
% around the controlled plant.
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also mask_Linear_Transfer_Function_DOB.

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

% Get linear dynamics matrices
[A, B, C, D, I, R, K] = jointObj.getDynamicsMatrices();

% Default values (off)
Q_td    = tf(0);
PQ_td	= tf(0);
PQ_ff	= tf(0);

% Check if DOB is enabled
if (DOB_enable_switch == 2)
    % If filename is given, assume the model and filters can be found in there
    if (loadData)
        if (isempty(modelFname))
            error('Error: Model/filters data file not found!');
        else
            % Load file and exit
            load(modelFname);
        end
    else
        % If no filename is given with model&filters from experimental
        % identification, do a linear simulation to obtain their
        % approximate values (this will significantly influence the results)

        [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = design_DOB_controller( jointObj, ...
                                                                Kp, Ki, Kd, N, 'parallel', ...
                                                                7, ... % Torque output
                                                                ff_comp_switch, ...
                                                                derivative_select, ...
                                                                f_c_FF, f_c_DOB, DOB_order  );
    end
else
    % If the DOB is disabled, make sure the prefilter is also disabled
    premult_enable_switch = 1;
end

% Discretize
% Using impulse-invariant transform for Q_td as that makes the filter not
% have feed-through, to avoid algebraic loops.
% This works well for simple low-pass filters. NOT advised for PQ_td or
% PQ_ff, as itresults in very wrong discretisation.
Q_td_d      = c2d(Q_td, jointObj.Ts, 'impulse');
PQ_td_d     = c2d(PQ_td, jointObj.Ts, 'tustin');
PQ_ff_d     = c2d(PQ_ff, jointObj.Ts, 'tustin');
