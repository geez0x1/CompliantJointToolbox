%% [A, B, C, I, D, K] = rigid(obj)
% Get dynamics matrices - fully rigid

function [A, B, C, I, D, K] = rigid(obj)
    
    % x = [q_b, q_b_dot]'

    % Inertia matrix
    I = obj.I_m + obj.I_g + obj.I_b;

    % Damping matrix
    d_m = obj.d_m;
    d_g = obj.d_g;
    d_b = obj.d_b; % shorthands %#ok<*PROP>
    D = d_m + d_g + d_b;

    % There is no stiffness (rigid joint)
    K = 0;

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];

    % Input
    k_t = obj.k_t;
    n   = obj.n;
    B	= [0, k_t*n/I(1,1)]';
    
        % Output
    C = [1, 0;  ... % motor position
         1, 0;  ... % gear position
         1, 0;  ... % link position
         0, 1;  ... % motor velocity
         0, 1;  ... % gear velocity
         0, 1;];... % link velocity
            
end
