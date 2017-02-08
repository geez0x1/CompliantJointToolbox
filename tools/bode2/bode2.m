% BODE2 Computes and plots magnitude and phase in frequency domain for a
% given transfer function. A wrapper around bode() and plot() that uses
% bode() for obtaining magnitude and phase info, and plot() for plotting. 
%
%   [ h, mag_db, phase, freq ] = bode2(P [, bodeOpt, opt, lineseries properties])
%
% Inputs::
%   P: transfer function to plot
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
%
% Examples::
%
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also bode_tuy, bode.

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

function [ h, mag_db, phase, freq ] = bode2( P, bodeOpt, opt, varargin )

    % Default arguments
    if (~exist('bodeOpt', 'var'))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var'))
        opt = bode2options;
    end

    % Obtain frequency, magnitude and phase data
    % First check whether XLim should be considered Hz or rad/s, based on
    % FreqUnits (bode() expects rad/s)
    if (strcmpi(bodeOpt.FreqUnits, 'Hz'))
        freq_range = {2*pi * bodeOpt.XLim{:}(1), 2*pi * bodeOpt.XLim{:}(2)};
    else
        freq_range = {bodeOpt.XLim{:}(1), 2*pi * bodeOpt.XLim{:}(2)};
    end
    [mag, phase, omega]	= bode(P, freq_range);
    mag = mag(:); phase = phase(:); omega = omega(:);
    mag_db	= mag2db(mag);
    f       = omega ./ (2*pi);
    
    % Set frequency units
    if (strcmpi(bodeOpt.FreqUnits, 'Hz'))
        freq = f;
    else
        freq = omega;
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

