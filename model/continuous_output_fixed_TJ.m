%% [A, B, C, I, D, K] = continuous_output_fixed_TJ(obj)
% Get dynamics matrices - output link fixed, torque-jerk states

function [A, B, C, I, D, K] = continuous_output_fixed_TJ(obj)
    
    % x = [tau_g, tau_l, tau_g_dot, tau_l_dot]'

    % Get position-velocity states
    [A, B, C, I, D, K] = obj.continuous_output_fixed();

    k_g = obj.k_g;
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    T = [   -k_g,	k_g;
            0,      -k_b    ];

    Tx = [	T,                  zeros(size(T));
            zeros(size(T)),  	T                   ];

    A = Tx*A/Tx;
    B = Tx*B;
    C = Tx\C;
    
end