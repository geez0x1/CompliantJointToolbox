% GETKALMAN Compute Kalman observer gains.
%
%   [kest, L, Cc] = getKalman(jointObj jOb, outputIdx,  var_u, var_y)
%
%  Calculate an optimal Kalman observer for the joint object jointObj with 
%  outputs specified by [outputIdx] and input and measurement variances var_u, var_y,
%  respectively. It returns the Kalman estimator kest, Kalman gain L, and
%  output matrix Cc.
%
% Inputs:
%   jointObj jOb: Joint object
%   outputIdx: Joint outputs measured by the Kalman filter
%   var_u: Input variance
%   var_y: Output variance
%
% Outputs:
%   kest: Kalman estimator dynamic system
%   L: Kalman gain
%   Cc: Output matrix that selects the outputs specified by outputIdx from
%       the output matrix C of the joint object.
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
% See also getObserver, getLinearDOB.

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

function [kest, L, Cc] = getKalman(jOb, outputIdx, var_u, var_y)
    %% Get state-space model
    sys     = jOb.getStateSpace();

    % Shorthands
    A       = sys.A;
    B       = sys.B;
    C       = sys.C;
    %D       = sys.D;

    % Create system with current input and outputs specified
    Ac   	= A;
    Bc   	= B(:,1);
    Cc    	= C(outputIdx,:);
    %Dc   	= zeros(size(Cc,1), 1);


    %% Design Kalman filter

    % x_dot	= Ax + Bu + Gw   	State equation
    % y  	= Cx + Du + Hw + v	Measurement equation

    % Build G, H
    G = Bc;                   	% Additive noise on the current (adding to u)
    H = zeros(size(Cc,1),1);	% No input noise feed-through

    % Construct sys_hat
    A_hat	= Ac;
    B_hat	= [Bc, G];
    C_hat	= Cc;
    D_hat	= [zeros(size(Cc,1),1), H];
    sys_hat = ss(A_hat, B_hat, C_hat, D_hat);

    % Define Kalman variance matrices
    Qn = var_u;                 % Input noise variance
    Rn = diag(var_y);           % Additive noise on the measurements

    % Calculate Kalman gains
    [kest, L, ~] = kalman(sys_hat, Qn, Rn);

    % Display results
    %disp('Calculated Kalman filter gain L:');
    %disp(L);
    %disp('and estimator state-space model kest.');
    %disp(kest);

end