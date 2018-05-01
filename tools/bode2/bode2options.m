% BODE2OPTIONS Default options for bode2
%
%   [ opt ] = bode2options()
%
% Inputs::
%
%
% Outputs::
%   opt: options struct
%
% Notes::
%
%
% Examples::
%
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also bode2.

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

function [ opt ] = bode2options()

    % General options
    opt                 = struct;
    opt.showTitle       = 1;
    opt.interpreter     = 'tex'; %#ok<*NASGU>
    opt.fontSize        = 8;
    opt.fig_w           = 600; % Default figure width [px]
    opt.fig_h           = 400; % Default figure height [px]
    opt.legendPos       = 'best'; % See doc on legend()
    opt.fixPhaseOffset  = 0; % Whether to fix the starting phase by removing multiples of 180 deg

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

