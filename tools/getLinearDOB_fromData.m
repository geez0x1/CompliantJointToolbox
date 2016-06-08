%% [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = getLinearDOB_fromData(jointName, t, u, y, f_c, id_Np, id_Nz, f_c_FF, f_c_DOB)
% This function calculates the approximated closed-loop transfer function
% Pc, low-pass Q-filters, and the inverted models for a DOB with premulti-
% plication control scheme. This is the experimental data version, that
% takes input-output data of the plant to approximate a plant model.
% This function returns the approximated closed-loop transfer function
% Pc, DOB filter Q_td, feed-forward filter Q_ff, DOB plant inversion +
% filter PQ_td, and feed-forward plant in version + filter PQ_ff.

function [Pc, Q_td, Q_ff, PQ_td, PQ_ff] = getLinearDOB_fromData(jointName, t, u, y, f_c, id_Np, id_Nz, f_c_FF, f_c_DOB)

    %% Default parameters
    if (~exist('f_c', 'var'))
        f_c = 60;    	% Model identification cut-off frequency [Hz]
    end
    if (~exist('id_Np', 'var'))
        id_Np = 4;       % Model number of poles []
    end
    if (~exist('id_Nz', 'var'))
        id_Nz = 1;       % Model number of zeros []
    end
    if (~exist('f_c_FF', 'var'))
        f_c_FF	= 40;	% Feed-forward cutoff frequency [Hz]
    end
    if (~exist('f_c_DOB', 'var'))
        f_c_DOB	= 60;	% DOB cutoff frequency [Hz]
    end

    % Bode options
    bodeOpt = bodeoptions;
    bodeOpt.FreqUnits = 'Hz';


    %% Get variables
    
    % Get joint object
    j = eval(jointName);
    
    % Cut-off frequencies
    omega_c_FF      = 2 * pi * f_c_FF;  % Feed-forward (model inv) LPF cutoff frequency [rad/s]
    omega_c_DOB     = 2 * pi * f_c_DOB; % DOB cutoff frequency [rad/s]
    %omega_c = 2 * pi * f_c;
    % The model identification cut-off is only used for displaying atm

    % Resample data to obtain uniform sampling for tfest()
    Ts      = j.Ts;         % Sampling time [s]
    t_RS    = 0:Ts:max(t);  % Resampled time
    u       = interp1(t, u, t_RS)';
    y       = interp1(t, y, t_RS)';
    t       = t_RS';

    % Plot bode plot of original data
    [f, mag_db, phase] = bode_tuy(t, u, y);

    
    %% Identification

    % Generate iddata object of data
    d = iddata(y, u, [], 'SamplingInstants', t);

    % Identify transfer function Pc
    Options = tfestOptions;
    Options.Display = 'on';
    Options.InitMethod = 'all';
    %Options.SearchMethod = 'gna';
    %Options.Focus = [0 30*2*pi];
    Pc = tfest(d, id_Np, id_Nz, Options);
    
    % Ensure Pc has unit DC gain
    %Pc = Pc / dcgain(Pc); % This is bad practice, but works if
    %dcgain(Pc)~1

    % Get magnitude and phase of Pc over f
    [mag_Pc, phase_Pc] = bode(Pc, 2*pi*f);
    mag_db_Pc	= mag2db(mag_Pc(:));
    phase_Pc	= phase_Pc(:);


    %% Plot original data and Pc approximation
    figure(1); clf;

    % Magnitude
    subplot(2,1,1); hold on;
    semilogx(f,mag_db);
    semilogx(f,mag_db_Pc, 'r');
    grid on
    xlim([0.1 f_c]);
    ylabel('Magnitude [dB]');
    legend('Experimental data', 'Pc');

    % Phase
    subplot(2,1,2); hold on;
    semilogx(f,phase);
    semilogx(f,phase_Pc, 'r');
    grid on;
    xlim([0.1 f_c]);
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');


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
    figure(2); clf; hold on;
    bode(Pc, bodeOpt);
    bode(inv(Pc), bodeOpt);
    bode(Q_td, bodeOpt);
    bode(Q_ff, bodeOpt);
    bode(PQ_td, bodeOpt);
    bode(PQ_ff, bodeOpt);
    xlim([0.1 100]);
    grid on;
    legend('P_c', 'P_c^{-1}', 'Q_{td}', 'Q_{ff}', 'PQ_{td}', 'PQ_{ff}');


    %% Save results
    saveName = [mfilename,'_results.mat'];
    save(saveName, 'Pc', 'Q_td', 'Q_ff', 'PQ_td', 'PQ_ff');
    disp(['Data saved to ', mfilename,'.mat']);
    
end
