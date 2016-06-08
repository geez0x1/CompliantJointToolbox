%% [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = Paine_analytical(jointName, Kp, Ki, Kd, N [, pid_form, ff_comp_switch, f_c_FF, f_c_DOB])
% This function calculates the approximated closed-loop transfer function
% Pc, low-pass Q-filters, and the inverted models for a DOB with premulti-
% plication control scheme. This is the standalone code, that does not
% require measurement data from a simulink model to design the model.
% Instead, it uses a linear state-space model to estimate the transfer
% function of the plant, thus the results are fundamentally different from
% the experimental version.
% This function returns the approximated closed-loop transfer function
% Pc, DOB filter Q_td, feed-forward filter Q_ff, DOB plant inversion +
% filter PQ_td, and feed-forward plant in version + filter PQ_ff.

function [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = Paine_analytical(jointName, Kp, Ki, Kd, N, pid_form, ff_comp_switch, f_c_FF, f_c_DOB)

    %% Default parameters
    if (~exist('pid_form', 'var'))
        pid_form = 'ideal';
    end
    if (~exist('ff_comp_switch', 'var'))
        ff_comp_switch = 1;     % Feed-forward/compensation (1=Compensation, 2=Feed-forward)
    end
    if (~exist('f_c_FF', 'var'))
        f_c_FF	= 40;           % Feed-forward cutoff frequency [Hz]
    end
    if (~exist('f_c_DOB', 'var'))
        f_c_DOB	= 60;           % DOB cutoff frequency [Hz]
    end

    % Bode options
    bodeOpt = bodeoptions;
    bodeOpt.FreqUnits = 'Hz';
    

    %% Get variables
    
    % Get joint object
    j = eval(jointName);

    
    % Control/system parameters
    k_b    	= j.k_b;                % Torsion bar stiffness [Nm/rad]
    n       = j.n;                  % Gearbox transmission ratio []
    k_t     = j.k_t;                % Torque constant [Nm/A]

    % Cut-off frequencies
    omega_c_FF  = 2 * pi * f_c_FF;  % Feed-forward (model inv) LPF cutoff frequency [rad/s]
    omega_c_DOB = 2 * pi * f_c_DOB; % DOB LPF cutoff frequency [rad/s]
    
    
    %% Get state-space system with torque output
    sys         = j.getStateSpace();
    sys         = ss(sys.A, sys.B, sys.C(2,:), 0);
    sys         = k_b * sys;  % Multiply by k_b to get torque output
    

    %% Build closed-loop system

    % PID controller
    if (strcmpi(pid_form, 'ideal'))
        Gf = pid(Kp, Kp*Ki, Kp*Kd, 1/N);	% ideal PID controller
    elseif (strcmpi(pid_form, 'parallel'))
        Gf = pid(Kp, Ki, Kd, 1/N);          % parallel PID controller
    else
        error('Invalid PID form');
    end
    Gf.u    = 'e';
    Gf.y    = 'pid_u';

    % See whether we need to build a closed-loop system with compensation or
    % feed-forward for the spring force
    if (ff_comp_switch == 1)        % Compensation
        % Compensation term
        C2 = 1 / (n * k_t);

        % Closed-loop transfer function
        P	= tf(sys);
        P2	= feedback(P, C2, +1); % Positive feedback
        Pf	= feedback(P2 * Gf, 1);

    elseif (ff_comp_switch == 2)    % Feed-forward
        % Feed-forward term
        C2      = tf(1 / (n * k_t));
        C2.u    = 'r';
        C2.y    = 'ff_u';

        % Closed-loop transfer function
        P       = tf(sys);
        P.u     = 'u';
        P.y     = 'y';
        sum_e   = sumblk('e = r - y');
        sum_ff  = sumblk('u = pid_u + ff_u');
        Pf      = connect(P, Gf, C2, sum_e, sum_ff, 'r', 'y');
        
    elseif (ff_comp_switch == 3)    % Spring force compensation off
        % Closed-loop transfer function
        P	= tf(sys);
        Pf	= feedback(P * Gf, 1);

    else
        error('Error: ff_comp_switch not set to 1, 2 or 3.');
        return;
    end
    
    
    %% Pc = Pf
    Pc = tf(Pf);


    %% Design low-pass Butterworth filters

    % Q_td
    [a, b] = butter(order(Pc), omega_c_DOB, 's');
    Q_td = tf(a,b);

    % Q_ff
    [a, b] = butter(order(Pc), omega_c_FF, 's');
    Q_ff = tf(a,b);

    % Pc^-1 * Q_td
    PQ_td = inv(Pc) * Q_td;

    % Pc^-1 * Q_ff
    PQ_ff = inv(Pc) * Q_ff;
    
    
    %% Show Bode plots of results
%     figure(2); clf; hold on;
%     bode(Gf,        bodeOpt);
%     bode(Pc,        bodeOpt);
%     bode(inv(Pc),	bodeOpt);
%     bode(Q_td,      bodeOpt);
%     bode(Q_ff,      bodeOpt);
%     bode(PQ_td,     bodeOpt);
%     bode(PQ_ff,     bodeOpt);
%     xlim([0.1 100]);
%     grid on;
%     legend('Gf', 'P_c', 'P_c^{-1}', 'Q_{td}', 'Q_{ff}', 'PQ_{td}', 'PQ_{ff}');


    %% Save results so we don't have to recalculate them all the time
    save('Paine_analytical_results.mat', 'Pc', 'Q_td', 'Q_ff', 'PQ_td', 'PQ_ff');
    disp('Data saved to Paine_analytical_results.mat');
    
end
