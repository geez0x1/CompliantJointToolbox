%% [f, mag_db, phase] = bode_tuyplot(t, u, y [,lineseries properties])
% Compute and plot magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.

function [f, mag_db, phase] = bode_tuyplot(t, u, y, varargin)
    
    % Get FFT
    [f, mag_db, phase] = bode_tuy(t, u, y);
    
    % Magnitude
    subplot(2,1,1); hold on;
    plot(f, mag_db,varargin{:});
    grid on
    ylabel('Magnitude [dB]');
    set(gca,'XScale','log');

    % Phase
    subplot(2,1,2); hold on;
    plot(f, phase,varargin{:});
    grid on;
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    set(gca,'XScale','log');
    
end