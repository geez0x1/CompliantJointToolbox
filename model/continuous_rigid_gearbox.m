%% [A, B, C, I, D, K] = continuous_rigid_gearbox(obj)
% Get dynamics matrices - rigid gearbox

function [A, B, C, I, D, K] = continuous_rigid_gearbox(obj)
    
    % x = [q_g, q_b, q_g_dot, q_b_dot]'

    % Inertia matrix
    I = diag([obj.I_m + obj.I_g, obj.I_b]);

    % Damping matrix
    d_m     = obj.d_m;
    d_g     = obj.d_g;
    d_b     = obj.d_b;
    d_gb	= obj.d_gb; % shorthands %#ok<*PROP>
    D = [   d_m + d_g + d_gb,       -d_gb; ...
            -d_gb,                  d_b + d_gb      ];

    % Stiffness matrix
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    K = [   k_b,        -k_b; ...
            -k_b,       k_b         ];

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];

    % Input
    k_t = obj.k_t;
    n   = obj.n;
    B	= [0, 0, k_t*n/I(1,1), 0]';
    
        % Output
    C = [1, 0, 0, 0;  ... % motor position
         1, 0, 0, 0;  ... % gear position
         0, 1, 0, 0;  ... % link position
         0, 0, 1, 0;  ... % motor velocity
         0, 0, 1, 0;  ... % gear velocity
         0, 0, 0, 1;];... % link velocity
            
end