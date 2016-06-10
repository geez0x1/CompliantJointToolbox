%% [A, B, C, I, D, K] = output_fixed_rigid_gearbox_TJ(obj)
% Get linear dynamics matrices - output link fixed, gearbox rigid, torque-jerk
% states

function [A, B, C, I, D, K] = output_fixed_rigid_gearbox_TJ(obj)
    
    % x = [tau_l, tau_l_dot]'

    % Get position-velocity states
    [A, B, C, I, D, K] = obj.output_fixed_rigid_gearbox();

    T = -obj.k_b;

    Tx = [	T,                  zeros(size(T));
            zeros(size(T)),     T               ];

    A = Tx*A/Tx;
    B = Tx*B;
    C = Tx\C;
    
end
