%% [A, B, C, I, D, K] = continuous_rigid_gearbox_TJ(obj)
% Get dynamics matrices - rigid gearbox - torque-jerk states

function [A, B, C, I, D, K] = continuous_rigid_gearbox_TJ(obj)
    
    % x = [tau_l, tau_e, tau_l_dot, tau_e_dot]'

    % Get position-velocity states
    [A, B, C, I, D, K] = obj.continuous_rigid_gearbox();

    k_b = obj.k_b;
    k_e = obj.k_e; % shorthands %#ok<*PROP>
    T = [   -k_b,	k_b;
            0,      k_e     ];

    Tx = [	T,                  zeros(size(T));
            zeros(size(T)),     T               ];

    A = Tx*A/Tx;
    B = Tx*B;
    C = Tx\C;
    
end