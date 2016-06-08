%% [kest, L, Cc] = getKalman(jointObj jOb, outputIdx,  var_u, var_y)
% Calculate an optimal Kalman observer for the jointObj jOb with outputs
% specified by [outputIdx]. It returns the Kalman estimator kest, Kalman
% gain L, and output matrix Cc. The input and measurement variances are
% given by var_u, var_y, respectively.

function [kest, L, Cc] = getKalman(jOb, outputIdx, var_u, var_y)
    %% Get state-space model
    sys     = jOb.getStateSpace();

    % Shorthands
    A       = sys.A;
    B       = sys.B;
    C       = sys.C;
    %D       = sys.D;

    % Create system with outputs specified
    Ac   	= A;
    Bc   	= B;
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