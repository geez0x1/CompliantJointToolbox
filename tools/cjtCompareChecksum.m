function flag = cjtCompareChecksum(input, ref)
%CJTCOMPARECHECKSUM Compare checksums of files or data.
%   This function is allows to compare checksums in one line of code
%   instead of repetedly compute checksums and perform the comparisons. It
%   uses cjtComputeChecksum to do that.
%
%   flag = cjtCompareChecksum(input, ref)
%
%   input:   File path, a numeric array or a character arry to be tested
%   ref:     Is a file path, numeric array or character array representing 
%            a precomputed checksum
%
%   'ref' 
%
%   flag:   True if checksums corresponding to 'input' and 'ref' are
%           matching.
%   
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also jointBuilder.

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

    chksum_1 = cjtComputeChecksum(input);

    cmp = 0;
    if ischar(ref)
        cmp = strcmp(chksum_1,ref); % ref can be a checksum itself, check if the checksum computed on the input matches
    end

    if ~cmp
        cmp = cjtComputeChecksum(input) == cjtComputeChecksum(ref); % compare element by element
        missmatches = find(cmp~= 1,1);
        flag = isempty(missmatches);  % does any element differ?
    else
        flag = cmp;
    end
end
