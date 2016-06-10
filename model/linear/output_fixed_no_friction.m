%% [A, B, C, I, D, K] = output_fixed(obj)
% Get dynamics matrices - output link fixed

function [A, B, C, I, D, K] = output_fixed(obj)
    
    % x = [q_m, q_g, q_m_dot, q_g_dot]'

    % Inertia matrix
    I = diag([obj.I_m , obj.I_g]);

    % Damping matrix
    d_m     = 0*obj.d_m;
    d_g     = 0*obj.d_g;
    d_mg    = obj.d_mg;
    d_gb	= obj.d_gb; % shorthands %#ok<*PROP>
    D = [	d_m + d_mg, 	-d_mg; ...
            -d_mg,      	d_g + d_mg + d_gb	];

    % Stiffness matrix
    k_g = obj.k_g;
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    K = [   k_g,    	-k_g; ...
            -k_g,    	k_g + k_b   ];

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];
        
    % Input
    k_t = obj.k_t;
    n	= obj.n;
    B	= [0, 0, k_t*n/I(1,1), 0]';
    
    % Output
    C = [1, 0, 0, 0;  ... % motor position
         0, 1, 0, 0;  ... % gear position
         0, 0, 0, 0;  ... % link position
         0, 0, 1, 0;  ... % motor velocity
         0, 0, 0, 1;  ... % gear velocity
         0, 0, 0, 0;];... % link velocity
    

end

