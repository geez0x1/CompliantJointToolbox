function [ opt ] = bodePlot2options()
%bodePlot2options Default options for bodePlot2

    opt                 = struct;
    opt.showTitle       = 1;
    opt.interpreter     = 'tex'; %#ok<*NASGU>
    opt.fontSize        = 8;
    opt.fig_w           = 600; % Default figure width [px]
    opt.fig_h           = 400; % Default figure height [px]
    opt.legendPos       = 'best'; % See doc on legend()

    % Default plot colours
    opt.C_R         = [1.0 0.4 0.4];
    opt.C_G         = [0.0 0.6 0.0];
    opt.C_B         = [0.3 0.3 1.0];
    opt.C_P         = [0.8 0.0 1.0];
    opt.C_R_light   = [1.0 0.6 0.6];
    opt.C_G_light   = [0.2 0.8 0.2];
    opt.C_B_light   = [0.6 0.6 1.0];
    opt.C_P_light   = [0.8 0.4 1.0];
    
    % Default colour order
    opt.C           = [opt.C_B; opt.C_R; opt.C_G; opt.C_P];
    opt.C_light     = [opt.C_B_light; opt.C_R_light; opt.C_G_light; opt.C_P_light];
    
end

