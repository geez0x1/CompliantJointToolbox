% BODE_TUYPLOT Compute and plot magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.
%
%   [freq, mag_db, phase] = bode_tuyplot(t, u, y [, resample, filter, bodeOpt, lineseries properties])
%
% Inputs::
%   t: time vector
%   u: input data vector
%   y: output data vector
%   resample: whether to resample data (for smaller filesize on export) (default: 0)
%   filter: whether to use zero-phase digital filtering for smoothing (default: 0)
%   bodeOpt: a bodeoptions struct (not all features are supported!)
%
% Outputs::
%   omega: frequency vector [rad/s]
%   mag_db: output magnitude vector in [db]
%   phase: output phase in [deg]
%
% Notes::
%   varargin holds additional plotting arguments passed to the plot
%   command.
%   The output values correspond to the plot: i.e. resampled and/or
%   filtered, and in frequency units as set in bodeOpt.
%
% Examples::
%
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also bode_tuy.

% Copyright (C) 2016, by Joern Malzahn, Wesley Roozing
%
% This file is part of the Compliant Joint Toolbox (CJT).
%
% CJT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CJT is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
% License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CJT. If not, see <http://www.gnu.org/licenses/>.
%
% For more information on the toolbox and contact to the authors visit
% <https://github.com/geez0x1/CompliantJointToolbox>

function [freq, mag_db, phase] = bode_tuyplot(t, u, y, resample, filter, bodeOpt, varargin)

    % Default arguments
    if (~exist('resample', 'var'))
        resample = 0;
    end
    if (~exist('filter', 'var'))
        filter = 0;
    end
    if (~exist('bodeOpt', 'var'))
        bodeOpt = bodeoptions;
    end
    
    % Get FFT
    [f, mag_db, phase] = bode_tuy(t, u, y);
    mag_db = mag_db(:); phase = phase(:); f = f(:);
    omega = f .* (2*pi);
    
    % Set frequency units
    if (strcmpi(bodeOpt.FreqUnits, 'Hz'))
        freq = f;
    else
        freq = omega;
    end
    
    % Resample data using logarithmic frequency space
    % Assume f starts at 0 (which it should always do)
    if (resample)
        a = log(freq(2)) / log(10);
        b = log(max(freq)) / log(10);
        n = round(max(freq) - min(freq));
        freq_RS   	= [0 logspace(a,b,n)];
        mag_db_RS  	= interp1(freq, mag_db, freq_RS);
        phase_RS 	= interp1(freq, phase, freq_RS);
        freq        = freq_RS;
        mag_db      = [mag_db_RS(1:end-1) mag_db(end)];
        phase       = [phase_RS(1:end-1) phase(end)];
    end
    
    % Filter magnitude and phase
    if (filter)
        windowSize	= 60;
        mag_db      = filtfilt(ones(1,windowSize) / windowSize, 1, mag_db);
        phase       = filtfilt(ones(1,windowSize) / windowSize, 1, phase);
    end
    
    % Depending on whether we show mag, phase, or both, do stuff
    if (strcmpi(bodeOpt.MagVisible, 'on') && strcmpi(bodeOpt.PhaseVisible, 'on'))
        % Magnitude
        subplot(2,1,1);
        semilogx(freq, mag_db, varargin{:});
        hold on;
        xlim(bodeOpt.XLim{1});
        ylabel('Magnitude [dB]', 'FontSize', bodeOpt.YLabel.FontSize, 'Interpreter', bodeOpt.YLabel.Interpreter);
        if (strcmp(bodeOpt.Grid, 'on'))
            grid on;
        end
        if (~strcmp(bodeOpt.Title.String, ''))
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end
        
        % Phase
        subplot(2,1,2);
        semilogx(freq, phase, varargin{:});
        hold on;
        xlim(bodeOpt.XLim{1});
        ylabel('Phase [deg]', 'FontSize', bodeOpt.YLabel.FontSize, 'Interpreter', bodeOpt.YLabel.Interpreter);
        if (strcmp(bodeOpt.FreqUnits, 'Hz'))
            xlabel('Frequency [Hz]', 'FontSize', bodeOpt.XLabel.FontSize, 'Interpreter', bodeOpt.XLabel.Interpreter);
        else
            xlabel('Frequency [rad/s]', 'FontSize', bodeOpt.XLabel.FontSize, 'Interpreter', bodeOpt.XLabel.Interpreter);
        end
        if (strcmp(bodeOpt.Grid, 'on'))
            grid on;
        end
        
    elseif (strcmpi(bodeOpt.MagVisible, 'on'))
        % Magnitude only
        semilogx(freq, mag_db, varargin{:});
        hold on;
        xlim(bodeOpt.XLim{1});
        ylabel('Magnitude [dB]', 'FontSize', bodeOpt.YLabel.FontSize, 'Interpreter', bodeOpt.YLabel.Interpreter);
        if (strcmp(bodeOpt.FreqUnits, 'Hz'))
            xlabel('Frequency [Hz]', 'FontSize', bodeOpt.XLabel.FontSize, 'Interpreter', bodeOpt.XLabel.Interpreter);
        else
            xlabel('Frequency [rad/s]', 'FontSize', bodeOpt.XLabel.FontSize, 'Interpreter', bodeOpt.XLabel.Interpreter);
        end
        if (strcmp(bodeOpt.Grid, 'on'))
            grid on;
        end
        if (~strcmp(bodeOpt.Title.String, ''))
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end
        
    elseif (strcmpi(bodeOpt.PhaseVisible, 'on'))
        % Phase only
        semilogx(freq, phase, varargin{:});
        hold on;
        xlim(bodeOpt.XLim{1});
        ylabel('Phase [deg]', 'FontSize', bodeOpt.YLabel.FontSize, 'Interpreter', bodeOpt.YLabel.Interpreter);
        if (strcmp(bodeOpt.FreqUnits, 'Hz'))
            xlabel('Frequency [Hz]', 'FontSize', bodeOpt.XLabel.FontSize, 'Interpreter', bodeOpt.XLabel.Interpreter);
        else
            xlabel('Frequency [rad/s]', 'FontSize', bodeOpt.XLabel.FontSize, 'Interpreter', bodeOpt.XLabel.Interpreter);
        end
        if (strcmp(bodeOpt.Grid, 'on'))
            grid on;
        end
        if (~strcmp(bodeOpt.Title.String, ''))
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end
        
    else
        error('Nothing to plot');
    end
    
end