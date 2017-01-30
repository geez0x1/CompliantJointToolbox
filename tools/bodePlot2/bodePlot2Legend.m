function [ h ] = bodePlot2Legend( data, bodeOpt, opt )
%bodePlot2Legend place legend in bodePlot2 Bode plots

    % Default arguments
    if (~exist('bodeOpt', 'var'))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var'))
        opt	= bodePlot2options;
    end
    
    % Switch to bottom subplot for plots with both mag and phase
    if (strcmp(bodeOpt.MagVisible, 'on') && strcmp(bodeOpt.PhaseVisible, 'on'))
        subplot(2,1,2);
    end
    
    % Place legend and set properties
    h = legend(data, 'Location', opt.legendPos);
    set(h, 'Interpreter', opt.interpreter);
    set(h, 'FontSize', opt.fontSize);

end

