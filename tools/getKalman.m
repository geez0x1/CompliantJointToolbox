% GETKALMAN Compute Kalman observer gains.
%
%   [kest, L, Cc] = getKalman(jointObj, inputIdx, outputIdx, var_u, var_y)
%
%  Calculate an optimal Kalman observer for the joint object jointObj with
%  inputs specified by [inputIdx], outputs specified by [outputIdx],
%  and input and measurement variances var_u, var_y, respectively. It
%  returns the Kalman estimator kest, Kalman gain L, and output matrix Cc.
%
% Inputs:
%   jointObj: Joint object
%   inputIdx: Joint inputs measured by the Kalman filter
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
%  For output_fixed models, when inputIdx does not contain the load
%  velocity (input 2), it is forcibly added, as it is required for the
%  system to be observable. Obviously during use this input can be zero for
%  truly locked actuator outputs, but the input needs to be given to
%  reconstruct the state of the model.
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

function [kest, L, Cc] = getKalman(jointObj, inputIdx, outputIdx, var_u, var_y)
    %% Get state-space model
    sys     = jointObj.getStateSpace();

    % Shorthands
    A       = sys.A;
    B       = sys.B;
    C       = sys.C;
    D       = sys.D;
    
    
    %% Check inputs
    
    % Output indices
    if (any(outputIdx < 1) || any(outputIdx > size(C,1)))
        error('getKalman error: Invalid output indices specified.');
    end
    
    % Output variance dimensions
    if (length(outputIdx) ~= length(var_y))
        error('getKalman error: Unequal number of outputs and output variances specified.');
    end
    
    % Input indices
    if (any(inputIdx < 1) || any(inputIdx > size(B,2)))
        error('getKalman error: Invalid input indices specified.');
    end
    
    % Add the second input (load motion) for output_fixed models
    % This is a bit of a hack which can possibly be removed with some
    % refactoring.
    if (~any(inputIdx == 2))
        if (    strcmp(jointObj.modelName, 'output_fixed')                              || ...
                strcmp(jointObj.modelName, 'output_fixed_no_friction')                  || ...
                strcmp(jointObj.modelName, 'output_fixed_rigid_gearbox')                || ...
                strcmp(jointObj.modelName, 'output_fixed_rigid_gearbox_no_friction')    )
            warning('getKalman warning: Adding secondary input (load velocity) as required input for Kalman filter for the system to be observable.');
            inputIdx = [inputIdx(:) 2];
        end
    end
    
    % Input variance dimensions
    if (length(inputIdx) ~= length(var_u))
        error('getKalman error: Invalid number of input variances specified.');
    end
    
    
    %% Create system with current input and outputs specified
    Ac      = A;
    Bc      = B(:,inputIdx);
    Cc      = C(outputIdx,:);
    Dc      = D(outputIdx,inputIdx);
    
    % Check observability
    % Using very small rank tolerance due to poles possibly being close to
    % the origin (effectively integrators)
    if (rank(obsv(Ac,Cc), eps) ~= length(A))
        error('getKalman error: This combination of model and measured outputs is not observable (observability matrix not full rank within tolerance of eps).');
    end


    %% Design Kalman filter

    % x_dot = Ax + Bu + Gw      State equation
    % y     = Cx + Du + Hw + v  Measurement equation

    % Build G, H
    G       = Bc;                           % Additive noise on the current (adding to u)
    H       = zeros(size(Cc,1),size(Bc,2)); % No input noise feed-through

    % Construct sys_hat
    A_hat   = Ac;
    B_hat   = [Bc, G];
    C_hat   = Cc;
    D_hat   = [Dc, H];
    sys_hat = ss(A_hat, B_hat, C_hat, D_hat);

    % Define Kalman variance matrices
    Qn      = diag(var_u);                  % Input noise variance
    Rn      = diag(var_y);                  % Additive noise on the measurements
    Nn      = zeros(length(var_u),length(var_y)); % Zero noise covariance

    % Check noise variances
    Z1      = [G zeros(size(G,1),size(H,1)); H eye(size(H,1))];
    Z2      = [Qn Nn;Nn' Rn];
    Z       = Z1 * Z2 * Z1';
    if (any(eig(Z)<0))
        error('getKalman error: Invalid noise variances (negative eigenvalues).');
    end
    
    % Calculate Kalman gains
    [kest, L, ~] = kalman(sys_hat, Qn, Rn);

    % Display results
    %disp('Calculated Kalman filter gain L:');
    %disp(L);
    %disp('and estimator state-space model kest.');
    %disp(kest);

end
