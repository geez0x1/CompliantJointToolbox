function [ outSum ] = cjtComputeChecksum( input)
%CJTCOMPUTECHECKSUM Compute checksum of files or data.
%   The implementation of this function reuses existing code. Since
%   the actual function used might change, this function is a wrapper to
%   decouple code.
%
%   outSum = cjtComputeChecksum(input)
%
%   input: numeric array, character array, file path
%   outSum: character array containing the checksum corresponding to the input
%
% Notes::
%   If 'input' is a path to an existing file on the path, the checksum will 
%   be computed for that file. Otherwise the checksum will be computed for
%   the corresponding character array.
%   
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

    if isnumeric(input)    % Numeric input?
        % Since we use Simlinki.getFileChecksum, we need to save the
        % numeric data to a file first. 
        fName = 'tmp_cjt_compute_checksum.txt';
        save(fName, 'input','-ascii');
        
        outSum = Simulink.getFileChecksum(fName);   
        
        % Delete temporary file.
        delete(fName);

    elseif ischar(input) && ~exist(input,'file')
        % The input is a character array that does not correspond to a file
        % name. Compute the checksum of the character array. 
        
        % Since we use Simlinki.getFileChecksum, we need to save the
        % numeric data to a file first. 
        fName = 'tmp_cjt_compute_checksum.txt';
        save(fName, 'input','-ascii');
        
        outSum = Simulink.getFileChecksum(fName);   
        
        % Delete temporary file.
        delete(fName);
        
    else
        % The input is a file path. We can use the built-in function directly. 
        outSum = Simulink.getFileChecksum(input);
    end

end
