% BODE_TUYPLOT Compute and plot magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.
%
%   [ h, mag_db, phase, freq ] = bode_tuyplot(t, u, y [, resample, filter, bodeOpt, opt, lineseries properties])
%
% Inputs::
%   t: time vector
%   u: input data vector
%   y: output data vector
%   resample: whether to resample data (for smaller filesize on export) (default: 0)
%   filter: whether to use zero-phase digital filtering for smoothing - implies resample=1 (default: 0)
%   bodeOpt: a bodeoptions struct (not all features are supported!)
%   opt: additional options (see bode2options)
%
% Outputs::
%   h: plot handle(s)
%   mag_db: output magnitude vector in [dB]
%   phase: output phase in [deg]
%   freq: frequency vector in chosen units
%
% Notes::
%   varargin holds additional plotting arguments passed to the plot
%   command.
%   The output values correspond to the plot: i.e. resampled or
%   resampled+filtered, and in frequency units as set in bodeOpt.
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

function [ h, mag_db, phase, freq ] = bode_tuyplot(t, u, y, resample, filter, bodeOpt, opt, varargin)

    % Default arguments
    if (~exist('resample', 'var') || isequal(resample,[]))
        resample = 0;
    end
    if (~exist('filter', 'var') || isequal(filter,[]))
        filter = 0;
    end
    if (~exist('bodeOpt', 'var') || isequal(bodeOpt,[]))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var') || isequal(opt,[]))
        opt = bode2options;
    end

    % If filtering is enabled, force resampling to be on as well
    % Otherwise the filtering window expressed in frequency is frequency
    % dependent.
    if (filter && ~resample)
        warning('[bode_tuyplot] Warning: Filtering enabled, forcing resampling.');
        resample = 1;
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

    % Treat the x-limits as the region of interest; that is, we remove the
    % other data (the frequency range of which is usually not sampled
    % anyway) so we can be faster and not influence the zero-phase
    % filtering below.
    idxs    = 1:length(freq);
    sIdx    = floor(interp1(freq, idxs, bodeOpt.XLim{1}(1)));
    eIdx    = ceil(interp1(freq, idxs, bodeOpt.XLim{1}(2)));
    freq    = freq(sIdx:eIdx);
    mag_db  = mag_db(sIdx:eIdx);
    phase   = phase(sIdx:eIdx);
    
    % Fix phase offset if requested by removing multiples of 180 deg at the
    % lowest frequency
    if (opt.fixPhaseOffset)
        phase = phase - round(phase(1)/180) * 180;
    end
    
    % Resample data using logarithmic frequency space.
    % freq may no longer start at 0 due to the ROI above.
    if (resample)
        a = log(freq(1)) / log(10);
        b = log(freq(end-1)) / log(10);
        n = 1000 * (b-a); % typically ~10x less data (for 2kHz sample time)
        freq_RS   	= logspace(a,b,n);
        mag_db_RS  	= interp1(freq, mag_db, freq_RS, 'pchip');
        phase_RS 	= interp1(freq, phase, freq_RS, 'pchip');
        freq        = freq_RS;
        mag_db      = [mag_db_RS(1:end-1) mag_db(end)];
        phase       = [phase_RS(1:end-1) phase(end)];
    end

    % Filter magnitude and phase
    if (filter)
        windowSize	= round(0.008 * n); % see n and generation of freq_RS
        mag_db      = filtfilt(ones(1,windowSize) / windowSize, 1, mag_db);
        phase       = filtfilt(ones(1,windowSize) / windowSize, 1, phase);
    end
    
    % Get figure handle and resize
    hFig = gcf;
    pos = get(hFig, 'Position');
    set(hFig, 'Position', [pos(1) pos(2) opt.fig_w opt.fig_h]);
    
    % Depending on whether we show mag, phase, or both, do stuff
    if (strcmpi(bodeOpt.MagVisible, 'on') && strcmpi(bodeOpt.PhaseVisible, 'on'))
        % Magnitude
        subplot(2,1,1);
        hMag = semilogx(freq, mag_db, varargin{:});
        hold on;
        xlim(bodeOpt.XLim{1});
        ylabel('Magnitude [dB]', 'FontSize', bodeOpt.YLabel.FontSize, 'Interpreter', bodeOpt.YLabel.Interpreter);
        if (strcmp(bodeOpt.Grid, 'on'))
            grid on;
        end
        if (opt.showTitle)
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end

        % Phase
        subplot(2,1,2);
        hPhase = semilogx(freq, phase, varargin{:});
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
        
        % Set handles
        h = [hMag, hPhase];
        
    elseif (strcmpi(bodeOpt.MagVisible, 'on'))
        % Magnitude only
        h = semilogx(freq, mag_db, varargin{:});
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
        if (opt.showTitle)
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end
        
    elseif (strcmpi(bodeOpt.PhaseVisible, 'on'))
        % Phase only
        h = semilogx(freq, phase, varargin{:});
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
        if (opt.showTitle)
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end
        
    else
        error('Nothing to plot');
    end

end

