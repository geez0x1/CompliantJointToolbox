% GETLINEARDOB Compute linear disturbance observer transfer functions.
%
%   [P, Q_td, PQ_td] = getLinearDOB(jointObj, omega_c, outputIdx [, DOB_order])
%
%   This function creates a DOB from a linear model jointObj with outputs
%   specified by [outputIdx]. It returns the plant model, Q-filter with
%   cut-off frequency omega_c, and the inverted plant + filter.
%
% Inputs::
%   jointObj: Joint object
%   omega_c: DOB Q-filter cut-off frequency in [rad/s]
%   outputIdx: Joint outputs measured by the observer
%   DOB_order: Order of the DOB (>= relative order of plant) (0 = auto)
%
% Outputs::
%   P: Plant transfer function
%   Q_td: Low-pass filter (Q-filter) with cut-off frequency omega_c
%   PQ_td: Inverted plant + filter
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
% See also getLinearDOB_fromData, getObserver, getKalman.

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

function [P, Q_td, PQ_td] = getLinearDOB(jointObj, omega_c, outputIdx, DOB_order)
    % Default parameters
    if (~exist('DOB_order', 'var') || isequal(DOB_order,[]))
        DOB_order = 0;
    end

    
    %% Get joint object and state-space system with current input and specified output
    sys     = jointObj.getStateSpace();
    sys     = ss(sys.A, sys.B(:,1), sys.C(outputIdx,:), sys.D(outputIdx,1));
    P       = tf(sys);

    % Check if P is a dynamic system
    % If we select the output position for a fixed-output plant for
    % example, P will be a static gain of 0, and a DOB cannot be
    % constructed.
    if (order(P) == 0)
        error('getLinearDOB error: Plant P has zero order for the chosen output, cannot continue.');
    end
    
    % Check selected DOB order
    % If it was set to zero before (no argument given), set it to the
    % relative order of P.
    if (DOB_order == 0)
        DOB_order = relativeOrder(P);
    end
    if (DOB_order < relativeOrder(P))
        error('DOB order needs to be equal to or larger than the relative order of the plant.');
    end
    
    %% Design low-pass Butterworth filters

    % Q_td
    [a, b]  = butter(DOB_order, omega_c, 's');
    Q_td    = tf(a,b);

    % P^-1 * Q_td
    PQ_td   = inv(P) * Q_td;
    
    
end
