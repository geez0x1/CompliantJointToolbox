%% [Pc, Q_td, PQ_td] = getLinearDOB(jointObj jOb, omega_c, outputIdx, doPlot)
% This function creates a DOB from a linear model jOb with outputs
% specified by [outputIdx]. It returns the plant model, Q-filter with
% cut-off frequency omega_c, and the inverted plant + filter. The doPlot
% flag plots Bode plots of the resulting transfer functions.

function [Pc, Q_td, PQ_td] = getLinearDOB(jOb, omega_c, outputIdx , doPlot)
    %% Get joint object and state-space system with 1 output
    sys     = jOb.getStateSpace();
    sys     = ss(sys.A, sys.B, sys.C(outputIdx,:), 0);
    Pc      = tf(sys);
    
	if ~exist('doPlot','var')
		doPlot = 0;
    end

    
    %% Design low-pass Butterworth filters

    % Q_td
    [a, b]	= butter(order(Pc), omega_c, 's');
    Q_td	= tf(a,b);

    % Pc^-1 * Q_td
    PQ_td	= inv(Pc) * Q_td;


    %% Show Bode plots of results
    % Bode options
    bodeOpt             = bodeoptions;
    bodeOpt.FreqUnits	= 'Hz';
    
    % Plot if required
	if doPlot
		figure(5); clf; hold on;
		bode(Pc, bodeOpt);
		bode(inv(Pc), bodeOpt);
		bode(Q_td, bodeOpt);
		bode(PQ_td, bodeOpt);
		xlim([0.1 100]);
		grid on;
		legend('P_c', 'P_c^{-1}', 'Q_{td}', 'PQ_{td}');
	end

    %% Save results so we don't have to recalculate them all the time
    save('DOB_results.mat', 'Pc', 'Q_td', 'PQ_td');
    disp('Saved Pc, Q_td, PQ_td to DOB_results.mat');
    
end
