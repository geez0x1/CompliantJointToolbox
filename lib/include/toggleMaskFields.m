%% toggleMaskFields( blockH, checkbox, fields [, onDirection ])
% Enable/disable mask fields based on the settings of checkboxes (supports
% multiple fields by using a cell for the 'fields' argument).
% Example: toggleMaskFields(gcbh, 'input_noise_enabled', 'var_u');
% Example: toggleMaskFields(gcbh, 'noise_disabled', {'var_u', 'var_y'}, 'off');

function [] = toggleMaskFields( blockH, checkbox, fields, onDirection )
    % Set default for onDirection
    if (~exist('onDirection', 'var'))
        onDirection     = 'on';     % By default, enabling the checkbox enables the fields
        offDirection	= 'off';
    elseif (~strcmpi(onDirection, 'on') && ~strcmpi(onDirection, 'off'))
        onDirection     = 'on';     % By default, enabling the checkbox enables the fields
        offDirection	= 'off';
    else
        if (strcmpi(onDirection, 'on'))
            offDirection = 'off';   % Checkbox off = fields off (default)
        else
            offDirection = 'on';    % Checkbox off = fields on
        end
    end

    % Get block mask names, mask enables
    MaskNames	= get(blockH, 'MaskNames');
    MaskEnables	= get(blockH, 'MaskEnables');

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

