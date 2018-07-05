function [outputArg] = OnOff2Bool(inputArg)
% ONOFF2BOOL Convert strings 'on' and 'off' to boolean equivalents.
% 
%   outputArg = OnOff2Bool(inputArg)
%
% Inputs::
%
%   inputArg:     String representing either 'on' or 'off'.
%
% Outputs::
%
%   outputArg: Boolean value 1 or 0.
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also Bool2OnOff.

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
switch lower(inputArg)
    case 'on'
        outputArg = boolean(1);
    case 'off'
        outputArg = boolean(0);
    otherwise
        error('Input argument must be ''on'' or ''off''!')
end

