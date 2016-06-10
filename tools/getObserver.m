% Calculate an observer for the joint object jOb with outputs specified by 
% [outputIdx]. The closed-loop poles based on LQR (with matrices Q, R) are
% multiplied by the place_gain to obtain the observer poles.
%
%   [sys_hat, L, Cc] = getObserver(jointObj jOb, outputIdx, place_gain, [Q, R])
%
% Inputs:
%   jointObj jOb: Joint object
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
% Notes::
%
%
% Examples::
%
%
% Author::
%  Joern Malzahn, jorn.malzahn@iit.it
%  Wesley Roozing, wesley.roozing@iit.it
%
% See also getKalman.

function [sys_hat, L, Cc] = getObserver(jOb, outputIdx, place_gain, Q, R)
    %% Parameters
    if ~exist('Q', 'var')
        Q = diag([0 1000 0 0]);
    end
    if ~exist('R', 'var')
        R = 1e-6;
    end

    
    %% Get state-space model
    sys     = jOb.getStateSpace();
    sys     = ss(sys.A, sys.B, eye(size(sys.A)), 0);

    % Shorthands
    A       = sys.A;
    B       = sys.B;
    C       = sys.C;
    D       = sys.D;

    % Create system with outputs specified
    Ac   	= A;
    Bc   	= B;
    Cc    	= C(outputIdx,:);
    Dc   	= 0;
    
    
    %% Check some dimensions
    if (size(Q) ~= size(Ac))
        error('getObserver error: size(Q) ~= size(Ac)');
        return;
    end
    if (length(R) ~= size(Bc,2))
        error('getObserver error: size(Q) ~= size(Ac)');
        return;
    end

    
    %% Design LQR controller

    % Calculate LQR gain matrix K for full state feedback
    [K_lqr, ~, ~] = lqr(sys, Q, R);

    % Calculate reference input premultiplication N
    % SISO, using position as output
    a           = [zeros(length(Bc),1); 1];
    N           = inv([Ac, Bc; Cc, Dc]) * a;
    N_x         = N(1:end-1);
    N_u         = N(end);
    N           = N_u + K_lqr * N_x;

    % Display results
    %disp('Calculated gain matrix K:');
    %disp(K_lqr);
    %disp('and premultiplication matrix N:');
    %disp(N);


    %% Design state observer

    % Obtain closed-loop controller poles for observer design
    cl_poles	= eig(A-B*K_lqr);
    obs_poles	= place_gain * min(real(cl_poles)) * ones(size(cl_poles));
    for i=1:length(obs_poles)
       obs_poles(i) = obs_poles(i) - place_gain*i; % Make poles distinct
    end

    % Calculate observer gain by placing poles
    L = place(Ac', Cc', obs_poles)';

    % Construct observer
    A_hat	= A;
    B_hat	= [B, L];
    C_hat	= eye(size(A));
    D_hat	= zeros(size(C_hat,1), size(B_hat,2));
    sys_hat = ss(A_hat, B_hat, C_hat, D_hat);
    
    % Display results
    %disp('Calculated observer gain L:');
    %disp(L);
    %disp('and observer state-space system sys_hat:');
    %sys_hat

end
