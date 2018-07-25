% GETOBSERVER Determine state observer gains
%
%   [sys_hat, L, Cc] = getObserver(jointObj, outputIdx, poles)
%
%  Calculate an observer for the joint object jointObj with outputs
%  specified by [outputIdx]. The closed-loop poles are selected by
%  the user as [poles].
%
% Inputs:
%   jointObj: Joint object
%   outputIdx: Joint outputs measured by the observer
%   place_gain: Pole placement gain
%   Q: LQR state weight
%   R: LQR input weight
%
% Outputs:
%   sys_hat: State observer dynamic system
%   L: Observer gain
%   Cc: Output matrix that selects the outputs specified by outputIdx from
%       the output matrix C of the joint object.
%
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

function [sys_hat, L, Cc] = getObserver(jointObj, outputIdx, poles)
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
        error('getObserver error: Invalid output indices specified.');
    end
    
    % Number of observer poles
    if (length(A) ~= length(poles))
        error(['getObserver error: Invalid number of observer poles specified (Expected ' num2str(length(A)) ').']);
    end
    
    
    %% Create system with inputs and outputs specified
    Ac      = A;
    Bc      = B(:,1);
    Cc      = C(outputIdx,:);
    %Dc      = D(outputIdx,1); % check this?
    
    % Check observability
    % Using very small rank tolerance due to poles possibly being close to
    % the origin (effectively integrators)
    if (rank(obsv(Ac,Cc), eps) ~= length(A))
        error('getObserver error: This combination of model and measured outputs is not observable (observability matrix not full rank within tolerance of eps).');
    end

    
    %% Design state observer
    
    % Check controllability for observer
    % Using very small rank tolerance due to poles possibly being close to
    % the origin (effectively integrators)
    if (rank(ctrb(Ac',Cc'), eps) ~= length(A))
        error('getObserver error: Observer matrices not controllable.');
    end

    % Calculate observer gain by placing poles
    L = place(Ac', Cc', poles)';
    
    % The commented-out code below designs the state observer based on LQR
%     %% Design LQR controller
%     
%     [K_lqr, ~] = getLQR(jointObj, outputIdx, Q, R);
% 
% 
%     %% Design state observer
% 
%     % Obtain closed-loop controller poles for observer design
%     cl_poles    = eig(A - Bc * K_lqr);
%     obs_poles   = place_gain * min(real(cl_poles)) * ones(size(cl_poles));
%     for i=1:length(obs_poles)
%        obs_poles(i) = obs_poles(i) * 1.1^(i-1); % Make poles distinct
%     end
% 
%     % Calculate observer gain by placing poles
%     L = place(Ac', Cc', obs_poles)';
% 
    % Construct observer
    A_hat   = A;
    B_hat   = [Bc, L];
    C_hat   = eye(size(A));
    D_hat   = zeros(size(C_hat,1), size(B_hat,2)); % check this?
    sys_hat = ss(A_hat, B_hat, C_hat, D_hat);
    

end
