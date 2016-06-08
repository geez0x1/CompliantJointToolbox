function [ tau ] = coulomb_asym(obj, x)
%COULOMB_ASYM Calculate asymmetric Coulomb friction torques
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
    
    % Build coefficient vector
    if (strcmp(obj.modelName, 'continuous_full'))
        c       = [0, 0, 0, obj.d_cm,	obj.d_cg,	obj.d_cb]';
        c_neg	= [0, 0, 0, obj.d_cm_n, obj.d_cg_n, obj.d_cb_n]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid_gearbox'))
        c       = [0, 0, obj.d_cm + obj.d_cg,       obj.d_cb]';
        c_neg   = [0, 0, obj.d_cm_n + obj.d_cg_n,	obj.d_cb_n]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed'))
        c       = [0, 0, obj.d_cm,      obj.d_cg]';
        c_neg   = [0, 0, obj.d_cm_n,    obj.d_cg_n]';
        
    elseif (strcmp(obj.modelName, 'continuous_output_fixed_rigid_gearbox'))
        c       = [0, obj.d_cm + obj.d_cg]';
        c_neg   = [0, obj.d_cm_n + obj.d_cg_n]';
        
    elseif (strcmp(obj.modelName, 'continuous_rigid'))
        c       = obj.d_cm + obj.d_cg + obj.d_cb;
        c_neg   = obj.d_cm_n + obj.d_cg_n + obj.d_cb_n;
        
    end
    
    % Build coefficients depending on signs of velocities
    cc(x>0) = c(x>0);
    cc(x<0) = c_neg(x<0);

    % If we have values over time per variable, pad the coefficients
    if (size(x,2) > 1)
        cc = padarray(cc, [1 size(x,2)-1], 'replicate', 'post');
    end
    
    % Calculate Coulomb friction torques
    tau = -cc .* atan(1000*x) / (pi/2);

end

