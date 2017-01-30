function [ ] = bodePlot2Title( bodeOpt, opt, paperMode, title_paperMode, title_normal )
%bodePlot2Title place title depending on paperMode and showTitle
    
    % Default arguments
    if (~exist('bodeOpt', 'var'))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var'))
        opt	= bodePlot2options;
    end
    
    % Check if we need to do stuff
    if (~opt.showTitle)
        return;
    end
    
    % Switch to top subplot for plots with both mag and phase
    if (strcmp(bodeOpt.MagVisible, 'on') && strcmp(bodeOpt.PhaseVisible, 'on'))
        subplot(2,1,1);
    end
    
    % Place title depending on paperMode
    if (paperMode)
        title(title_paperMode, 'Interpreter', bodeOpt.Title.Interpreter);
    else
        title(title_normal, 'Interpreter', bodeOpt.Title.Interpreter);
    end
end

