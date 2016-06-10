% Compute and plot magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.
%
%   [f, mag_db, phase] = bode_tuyplot(t, u, y [,lineseries properties])
%
% Inputs:
%   t: time vector
%   u: input data vector
%   y: output data vector
%
% Outputs:
%   f: frequency vector
%   mag_db: output magnitude vector in [db]
%   phase: output phase in [deg]
%
% Notes::
%   varargin holds additional plotting arguments passed to the plot command
%
% Examples::
%
%
% Author::
%  Joern Malzahn, jorn.malzahn@iit.it
%  Wesley Roozing, wesley.roozing@iit.it
%
% See also bode_tuy.

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