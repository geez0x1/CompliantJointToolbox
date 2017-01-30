function [ ] = resizeFig( h, width, height )
%resizeFig Resize the figure specified by the handle h to the width and
% height specified, without moving it.

    pos     = get(h, 'Position');
    xPos	= pos(1);
    yPos	= pos(2);
    set(h, 'Position', [xPos yPos width height]);

end

