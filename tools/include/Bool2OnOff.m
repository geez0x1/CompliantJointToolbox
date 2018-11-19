function [outputArg] = Bool2OnOff(inputArg)
% ONOFF2BOOL Convert booleans to equivalent strings 'on' and 'off'.
% 
%   outputArg = Bool2OnOff(inputArg)
%
% Inputs::
%
%   inputArg:    Boolean value 1 or 0.
%
% Outputs::
%
%   outputArg:  String representing either 'on' or 'off'.
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also OnOff2Bool.

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

switch (inputArg)
    case 0
        outputArg = 'off';
    case 1
        outputArg = 'on';
    otherwise
        error('Input argument must be either 1 or 0.')
end