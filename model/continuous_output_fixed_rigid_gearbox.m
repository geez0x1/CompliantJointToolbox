%% [A, B, C, I, D, K] = continuous_output_fixed_rigid_gearbox(obj)
% Get dynamics matrices - output link fixed, gearbox rigid

function [A, B, C, I, D, K] = continuous_output_fixed_rigid_gearbox(obj)
    
    % x = [q_g, q_g_dot]'

    % Inertia matrix
    I = obj.I_m + obj.I_g;

    % Damping matrix
    D = obj.d_m + obj.d_g + obj.d_gb;

    % Stiffness matrix
    K = obj.k_b;

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
         0, 0;  ... % link position
         0, 1;  ... % motor velocity
         0, 1;  ... % gear velocity
         0, 0;];... % link velocity

end

