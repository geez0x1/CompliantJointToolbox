% BODE2LEGEND Place legend in bode2 Bode plots
%
%   [ h ] = bode2legend( data, bodeOpt, opt )
%
% Inputs::
%   data: legend strings [cell]
%   bodeOpt: a bodeoptions struct (not all features are supported!)
%   opt: additional options (see bode2options)
%
% Outputs::
%   h: legend handle
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

function [ h ] = bode2legend( data, bodeOpt, opt )

    % Default arguments
    if (~exist('bodeOpt', 'var') || isequal(bodeOpt,[]))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var') || isequal(opt,[]))
        opt	= bode2options;
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

