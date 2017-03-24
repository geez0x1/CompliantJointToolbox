% GETOPENLOOPDOBTF Compute the transfer function resulting from
% interconnecting a plant P with DOB enforcing a nominal plant model Pn.
%
%   [H_dob, Q_td, PQ_td] = getOpenLoopDOBTF(jointObj_P, jointObj_Pn, omega_c, outputIdx, doPlot)
%
%   This function creates a DOB from the nominal plant model jointObj_Pn,
%   and closes it around the 'real' plant model P, based on the output
%   specified by outputIdx. It returns the input-output transfer function
%   of the plant + DOB, the DOBs Q-filter with cut-off frequency omega_c,
%   and the inverted plant + filter used in the DOB. The doPlot flag plots
%   Bode plots of the resulting transfer functions.
%
% Inputs::
%   jointObj_P: 'real' plant joint object
%   jointObj_Pn: Nominal plant joint object
%   omega_c: DOB Q-filter cut-off frequency in [rad/s]
%   outputIdx: Joint outputs measured by the observer
%   doPlot: Whether to plot the results
%
% Outputs::
%   H_dob: Input-output transfer function of the plant+DOB
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
% See also getLinearDOB, getLinearDOB_fromData, getObserver, getKalman.

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

function [H_dob, Q_td, PQ_td] = getOpenLoopDOBTF(jointObj_P, jointObj_Pn, omega_c, outputIdx, doPlot)
    % Default parameters
    if (~exist('doPlot','var') || isequal(doPlot,[]))
        doPlot = 0;
    end

    
    %% Get joint objects and state-space systems with current input and specified output
    sys_P   = jointObj_P.getStateSpace();
    sys_P   = ss(sys_P.A, sys_P.B(:,1), sys_P.C(outputIdx,:), sys_P.D(outputIdx,1));
    P       = tf(sys_P);
    
    sys_Pn  = jointObj_Pn.getStateSpace();
    sys_Pn	= ss(sys_Pn.A, sys_Pn.B(:,1), sys_Pn.C(outputIdx,:), sys_Pn.D(outputIdx,1));
    Pn      = tf(sys_Pn);

    % Check if P and Pn are dynamic systems
    % If we select the output position for a fixed-output plant for
    % example, Pn will be a static gain of 0, and a DOB cannot be
    % constructed.
    if (order(P) == 0 || order(Pn) == 0)
        error('getOpenLoopDOBTF error: Either P or Pn has zero order for the chosen output, cannot continue.');
    end
    
    
    %% Design low-pass Butterworth filter and nominal plant inversion

    [Q_td, PQ_td] = getLinearDOB(jointObj_Pn, omega_c, outputIdx, 0); % don't plot here, as we also plot below
    
    
    %% Calculate full loop of plant + DOB

    H_dob = P / (1 + Q_td * (P/Pn - 1));
    

    %% Show Bode plots of results
    % Bode options
    bodeOpt             = bodeoptions;
    bodeOpt.FreqUnits   = 'Hz';
    
    % Plot if required
    if doPlot
        figure(5); clf; hold on;
        bode(P, bodeOpt);
        bode(Pn, bodeOpt);
        bode(H_dob);
        bode(inv(Pn), bodeOpt);
        bode(Q_td, bodeOpt);
        bode(PQ_td, bodeOpt);
        xlim([0.1 100]);
        grid on;
        legend('P', 'P_n', 'H_{dob}', 'P_n^{-1}', 'Q_{td}', 'PQ_{td}');
    end
    
end
