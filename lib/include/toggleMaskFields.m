% TOGGLEMASKFIELDS Enable/disable mask fields based on the settings of checkboxes
% (supports multiple fields by using a cell for the 'fields' argument).
%
% toggleMaskFields( blockH, checkbox, fields [, onDirection ])
%
% Inputs::
%
% blockH:      handle to the Simulink block under consideration
% checkbox:    handle to the toggled checkbox
% fields:      cell array of fields in the mask 
% onDirection: Either 'on' or 'off'. By default, enabling the checkbox
%              enables the fields ('on')
%
% Examples::
%  1)
%   toggleMaskFields(gcbh, 'input_noise_enabled', 'var_u');
%  2)
%   toggleMaskFields(gcbh, 'noise_disabled', {'var_u', 'var_y'}, 'off');
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also mask_LQR_Full_Measured_State.

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

function [] = toggleMaskFields( blockH, checkbox, fields, onDirection )
    % Set default for onDirection
    if (~exist('onDirection', 'var'))
        onDirection     = 'on';     % By default, enabling the checkbox enables the fields
        offDirection    = 'off';
    elseif (~strcmpi(onDirection, 'on') && ~strcmpi(onDirection, 'off'))
        onDirection     = 'on';     % By default, enabling the checkbox enables the fields
        offDirection    = 'off';
    else
        if (strcmpi(onDirection, 'on'))
            offDirection = 'off';   % Checkbox off = fields off (default)
        else
            offDirection = 'on';    % Checkbox off = fields on
        end
    end

    % Get block mask names, mask enables
    MaskNames   = get(blockH, 'MaskNames');
    MaskEnables = get(blockH, 'MaskEnables');

    % Find the index/indices of the field(s)
    if (iscell(fields)) % Multiple fields
        j = 1; % Valid fields counter
        for i=1:numel(fields)
            fieldIdx = find(strcmpi(MaskNames, fields(i)));
            if (numel(fieldIdx)>0)
                idx(j) = fieldIdx;
                j = j+1;
            end
        end
    else % Single field
        fieldIdx = find(strcmpi(MaskNames, fields));
        if (numel(fieldIdx)>0)
            idx = fieldIdx;
        else
            error(['toggleMaskFields error: Unable to find field ' fields]);
        end
    end

    % If the checkbox is enabled, enable the fields
    if (strcmpi(get(blockH, checkbox), 'on'))
        for i=1:numel(idx)
            MaskEnables{idx(i)} = onDirection;
        end
    else
        for i=1:numel(idx)
            MaskEnables{idx(i)} = offDirection;
        end
    end
    
    % Actually set MaskEnables field into block
    set(blockH, 'MaskEnables', MaskEnables);
end

