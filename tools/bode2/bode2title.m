% BODE2TITLE Place title depending on paperMode and showTitle
%
%   [ ] = bode2title( bodeOpt, opt, paperMode, title_paperMode, title_normal )
%
% Inputs::
%   bodeOpt: a bodeoptions struct (not all features are supported!)
%   opt: additional options (see bode2options)
%   paperMode: flag that determines which title to use
%   title_paperMode: paper mode title to use
%   title_normal: normal title to use (optional, otherwise equal to
%                 title_paperMode)
%
% Outputs::
%
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

function [ ] = bode2title( bodeOpt, opt, paperMode, title_paperMode, title_normal )
    
    % Default arguments
    if (~exist('bodeOpt', 'var'))
        bodeOpt = bodeoptions;
    end
    if (~exist('opt', 'var'))
        opt	= bode2options;
    end
    if (~exist('title_normal', 'var'))
        title_normal = title_paperMode;
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

