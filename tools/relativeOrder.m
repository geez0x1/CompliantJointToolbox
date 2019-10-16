% RELATIVEORDER Compute relative order of transfer function
% =========================================================================
%
% [relOrder] = relativeOrder(P)
%
%  Description:
%    Time and signal vector must have equal length. Time steps are assumed
%    to be equidistant.
%
%  Inputs::
%    P: Transfer function
%
%  Outputs::
%    relOder: Relative order of P []
%
%
%  Authors:
%    Joern Malzahn
%    Wesley Roozing
%
%  See also fft, ifft.
%
% Copyright (C) 2019, by Joern Malzahn, Wesley Roozing
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
%

function [ relOrder ] = relativeOrder( P )
    % Check whether P is actually a transfer function
    if (~isequal(P, tf(P)))
        error('Error: P is not a transfer function');
    end

    % Get numerator and denominator orders
    [~, a] = find(P.num{:} ~= 0, 1, 'first');
    [~, b] = find(P.den{:} ~= 0, 1, 'first');
    numOrder = length(P.num{:}) - a;
    denOrder = length(P.den{:}) - b;
    
    % Note we could also use the padding in front of the numerator
    % coefficient array (as MATLAB keeps them equal length), but since this
    % is a bit of a quirk we'll use this neater method.
    
    % Compute relative order
    relOrder = denOrder - numOrder;

end

