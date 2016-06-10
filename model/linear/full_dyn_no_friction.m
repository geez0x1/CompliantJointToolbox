%% [A, B, C, I, D, K] = full_dyn_no_friction(obj)
% Get dynamics matrices - default, no friction

function [A, B, C, I, D, K] = full_dyn_no_friction(obj)
    
    % x = [q_m, q_g, q_b, q_m_dot, q_g_dot, q_b_dot]'

    % Inertia matrix
    I = diag([obj.I_m, obj.I_g, obj.I_b]);

    % Damping matrix
    d_m     = 0;
    d_g     = 0*obj.d_g;
    d_b     = 0*obj.d_b;
    d_mg    = obj.d_mg;
    d_gb    = obj.d_gb; % shorthands %#ok<*PROP>
    D = [   d_m + d_mg,     -d_mg,                  0; ...
            -d_mg,          d_g + d_mg + d_gb,      -d_gb; ...
            0,              -d_gb,                  d_b + d_gb	];

    % Stiffness matrix
    k_g = obj.k_g;
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    K = [   k_g,            -k_g,           0; ...
            -k_g,           k_g + k_b,      -k_b; ...
            0,              -k_b,           k_b         ];

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];

    % Input
    k_t = obj.k_t;
    n   = obj.n;
    B	= [0, 0, 0, k_t*n/I(1,1), 0, 0]';
    
    % Output
    C = eye(size(A,2));

end
