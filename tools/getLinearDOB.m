% GETLINEARDOB Compute linear disturbance observer transfer functions.
%
%   [Pc, Q_td, PQ_td] = getLinearDOB(jointObj jOb, omega_c, outputIdx, doPlot)
%
%   This function creates a DOB from a linear model jOb with outputs
%   specified by [outputIdx]. It returns the plant model, Q-filter with
%   cut-off frequency omega_c, and the inverted plant + filter. The doPlot
%   flag plots Bode plots of the resulting transfer functions.
%
% Inputs::
%   jointObj: Joint object
%   omega_c: DOB Q-filter cut-off frequency in [rad/s]
%   outputIdx: Joint outputs measured by the observer
%   doPlot: Whether to plot the results
%
% Outputs::
%   Pc: Plant transfer function
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

function [Pc, Q_td, PQ_td] = getLinearDOB(jOb, omega_c, outputIdx , doPlot)
    %% Get joint object and state-space system with 1 output
    sys     = jOb.getStateSpace();
    sys     = ss(sys.A, sys.B, sys.C(outputIdx,:), 0);
    Pc      = tf(sys);
    
	if ~exist('doPlot','var')
		doPlot = 0;
    end

    
    %% Design low-pass Butterworth filters

    % Q_td
    [a, b]	= butter(order(Pc), omega_c, 's');
    Q_td	= tf(a,b);

    % Pc^-1 * Q_td
    PQ_td	= inv(Pc) * Q_td;


    %% Show Bode plots of results
    % Bode options
    bodeOpt             = bodeoptions;
    bodeOpt.FreqUnits	= 'Hz';
    
    % Plot if required
	if doPlot
		figure(5); clf; hold on;
		bode(Pc, bodeOpt);
		bode(inv(Pc), bodeOpt);
		bode(Q_td, bodeOpt);
		bode(PQ_td, bodeOpt);
		xlim([0.1 100]);
		grid on;
		legend('P_c', 'P_c^{-1}', 'Q_{td}', 'PQ_{td}');
	end

    %% Save results so we don't have to recalculate them all the time
    save('DOB_results.mat', 'Pc', 'Q_td', 'PQ_td');
    disp('Saved Pc, Q_td, PQ_td to DOB_results.mat');
    
end
