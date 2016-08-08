%__________________________________________________________________
% Command-line yes/no prompts
function [out] = confirm(prompt, default)
    ansr = 'temp';

    % Keep asking the prompt until we get '', 'n', or 'y'
    while (~strcmpi(ansr, 'y') && ~strcmpi(ansr, 'n') && ~strcmp(ansr, ''))
        ansr = input(prompt, 's');
    end

    % Check answer. If default==1 then accept '' as yes
    if ((default == 1 && strcmp(ansr, '')) || strcmpi(ansr, 'y'))
        out = 1;
        return;
    end

    % No otherwise
    out = 0;
end