function [ tau ] = coulomb(obj, x)
%COULOMB Calculate Coulomb friction torques
% Arguments:
% jointObj obj
% state x

    % x = [q_m; q_g; q_b; q_m_dot; q_g_dot, q_b_dot'];  continuous_full
    % x = [q_g, q_b, q_g_dot, q_b_dot]'                 continuous_rigid_gearbox
    % x = [q_m, q_g, q_m_dot, q_g_dot]'                 continuous_output_fixed
    % x = [q_g, q_g_dot]'                               continuous_output_fixed_rigid_gearbox
    
    % Preallocate coefficient vector
    c = zeros(size(x));
    
    % Build coefficient vector
    if (strcmp(obj.modelName, 'continuous_full'))
        c = [0, 0, 0, obj.d_cm, obj.d_cg, obj.d_cb]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid_gearbox'))
        c = [0, 0, obj.d_cm + obj.d_cg, obj.d_cb]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed'))
        c = [0, 0, obj.d_cm, obj.d_cg]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed_rigid_gearbox'))
        c = [0, obj.d_cm + obj.d_cg]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid'))
        c = obj.d_cm + obj.d_cg + obj.d_cb;
        
    end

    % If we have values over time per variable, pad the coefficients
    if (size(x,2) > 1)
        c = padarray(c, [1 size(x,2)-1], 'replicate', 'post');
    end

    % Calculate Coulomb friction torques
    tau = -c .* atan(1000*x) / (pi/2);

end

