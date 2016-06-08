function [ tau ] = viscous_asym(obj, x)
%VISCOUS_ASYM Calculate asymmetric viscous friction torque
% Arguments:
% jointObj obj
% state x

    % x = [q_m; q_g; q_b; q_m_dot; q_g_dot, q_b_dot'];  continuous_full
    % x = [q_g, q_b, q_g_dot, q_b_dot]'                 continuous_rigid_gearbox
    % x = [q_m, q_g, q_m_dot, q_g_dot]'                 continuous_output_fixed
    % x = [q_g, q_g_dot]'                               continuous_output_fixed_rigid_gearbox
    
    % Preallocate coefficient vectors
    c       = zeros(size(x));
    c_neg	= zeros(size(x));
    cc      = zeros(size(x));
    
    % Get some shorthands
    d_m     = obj.d_m;
    d_g     = obj.d_g;
    d_b     = obj.d_b;
    d_m_n   = obj.d_m_n;
    d_g_n   = obj.d_g_n;
    d_b_n   = obj.d_b_n;
    
    % Build coefficient vector
    % Because the positive part is built into the linear dynamics, we
    % compensate for them in this asymmetric nonlinear model to obtain the
    % desired damping values. For the same reason, the positive part is zero.
    if (strcmp(obj.modelName, 'continuous_full'))
        c       = [0, 0, 0, 0,              0,              0               ]';
        c_neg	= [0, 0, 0, -d_m + d_m_n,	-d_g + d_g_n,	-d_b + d_b_n	]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid_gearbox'))
        c       = [0, 0, 0,                             0               ]';
        c_neg   = [0, 0, -d_m - d_g + d_m_n + d_g_n,	-d_g + d_g_n	]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed'))
        c       = [0, 0, 0,             0               ]';
        c_neg   = [0, 0, -d_m + d_m_n,	-d_g + d_g_n	]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed_rigid_gearbox'))
        c       = [0, 0                             ]';
        c_neg   = [0, -d_m - d_g + d_m_n + d_g_n	]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid'))
        c       = 0;
        c_neg   = -d_m - d_g - d_b + d_m_n + d_g_n + d_b_n;
        
    end
    
    % Build coefficients depending on signs of velocities
    cc(x>0) = c(x>0);
    cc(x<0) = c_neg(x<0);

    % If we have values over time per variable, pad the coefficients
    if (size(x,2) > 1)
        cc = padarray(cc, [1 size(x,2)-1], 'replicate', 'post');
    end
    
    % Calculate asymmetric viscous friction torques
    tau = -cc .* x;

end

