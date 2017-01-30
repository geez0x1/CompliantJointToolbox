function [ h, mag_db, phase, freq ] = bodePlot2( P, bodeOpt, opt, varargin )
%bodePlot2 A wrapper around bode() and plot() that uses bode() for
%obtaining magnitude and phase info, and plot() for plotting. Because
%MATLAB and its stupid Bode() function suck.

    % Default arguments
    if (~exist('bodeOpt', 'var'))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var'))
        opt = bodePlot2options;
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
    if (strcmp(bodeOpt.FreqUnits, 'Hz'))
        freq = f;
    else
        freq = omega;
    end
    
    % Get figure handle and resize
    h = gcf;
    pos = get(h, 'Position');
    set(h, 'Position', [pos(1) pos(2) opt.fig_w opt.fig_h]);
    
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
        if (opt.showTitle)
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
        if (opt.showTitle)
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
        if (opt.showTitle)
            title(bodeOpt.Title.String, 'FontSize', bodeOpt.Title.FontSize, 'Interpreter', bodeOpt.Title.Interpreter);
        end
        
    else
        error('Nothing to plot');
    end

end

