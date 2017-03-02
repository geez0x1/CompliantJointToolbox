% GETOBSERVER Determine state observer gains
%
%   [sys_hat, L, Cc] = getObserver(jointObj, outputIdx, place_gain, [Q, R])
%
%  Calculate an observer for the joint object jointObj with outputs specified by 
%  [outputIdx]. The closed-loop poles based on LQR (with matrices Q, R) are
%  multiplied by the place_gain to obtain the observer poles.
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

function [sys_hat, L, Cc] = getObserver(jointObj, outputIdx, place_gain, Q, R)
    %% Parameters
    if (~exist('Q', 'var') || isequal(Q,[]))
        Q = diag([0 1000 0 0]);
    end
    if (~exist('R', 'var') || isequal(R,[]))
        R = 1e-6;
    end
    if (length(outputIdx) > 1)
        error('getLQR error: More than one output weight specified.');
    end

    
    %% Get state-space model
    sys     = jointObj.getStateSpace();
    sys     = ss(sys.A, sys.B, sys.C, sys.D);

    % Shorthands
    A       = sys.A;
    B       = sys.B(:,1); % Use only the current input
    C       = sys.C;
    D       = sys.D(:,1); % Use only the current input

    % Create system with current input and outputs specified
    Ac      = A;
    Bc      = B;
    Cc      = C(outputIdx,:);
    Dc      = D(outputIdx,1);
    sys     = ss(Ac, Bc, Cc, Dc);
    
    
    %% Check some dimensions
    if (size(Q) ~= size(Ac))
        error('getObserver error: size(Q) ~= size(Ac)');
    end
    if (length(R) ~= size(Bc,2))
        error('getObserver error: size(R) ~= size(Bc)');
    end

    
    %% Design LQR controller
    
    [K_lqr, ~] = getLQR(jointObj, outputIdx, Q, R);


    %% Design state observer

    % Obtain closed-loop controller poles for observer design
    cl_poles    = eig(A - Bc * K_lqr);
    obs_poles   = place_gain * min(real(cl_poles)) * ones(size(cl_poles));
    for i=1:length(obs_poles)
       obs_poles(i) = obs_poles(i) * 1.1^(i-1); % Make poles distinct
    end

    % Calculate observer gain by placing poles
    L = place(Ac', Cc', obs_poles)';

    % Construct observer
    A_hat   = A;
    B_hat   = [Bc, L];
    C_hat   = eye(size(A));
    D_hat   = zeros(size(C_hat,1), size(B_hat,2));
    sys_hat = ss(A_hat, B_hat, C_hat, D_hat);
    
    % Display results
    %disp('Calculated observer gain L:');
    %disp(L);
    %disp('and observer state-space system sys_hat:');
    %sys_hat

end
