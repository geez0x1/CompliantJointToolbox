function [ tau ] = viscous_pos_dep(obj, x)
%VISCOUS_POS_DEP Calculate position-dependent viscous friction
% Arguments:
% jointObj obj
% state x

    % x = [q_m; q_g; q_b; q_m_dot; q_g_dot, q_b_dot'];  continuous_full
    % x = [q_g, q_b, q_g_dot, q_b_dot]'                 continuous_rigid_gearbox
    % x = [q_m, q_g, q_m_dot, q_g_dot]'                 continuous_output_fixed
    % x = [q_g, q_g_dot]'                               continuous_output_fixed_rigid_gearbox
    
    % THIS FILE IS A PROOF OF CONCEPT AND NEEDS TO BE FINISHED
    
    % Temp
    a = 1e-6;   % Amplitude
    d = 1/2;    % Frequency
    p = pi;     % Phase
    
    c = zeros(size(x));
    
    % Build coefficient vector
    if (strcmp(obj.modelName, 'continuous_full'))
        c = [0, 0, 0, ...
                a * sin(x(1)/d + p), ...
                a * sin(x(2)/d + p), ...
                a * sin(x(3)/d + p)       ]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid_gearbox'))
        c = [0, 0, ...
                a * sin(x(1)*d + p) + a * sin(obj.n*x(1)*d + p) * obj.n, ...
                a * sin(x(2)*d + p)       ]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed'))
        c = [0, 0, ...
                a * sin(x(1)/d + p), ...
                a * sin(x(2)/d + p)       ]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed_rigid_gearbox'))
        c = [0, a * sin(x(1)*d + p) + a * sin(obj.n*x(1)*d + p) * obj.n]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid'))
        c = a * sin(x(1)*d + p) + a * sin(obj.n*x(1)*d + p) * obj.n;
        
    end

    % Calculate position-dependent viscous friction torques
    tau = -c .* x;

end

