function printpdf(figHandle,outfilename,varargin)
% PRINTPDF A wrapper to the built-in function 'print'.
%
%   printpdf(figHandle,fileName[,options])
%
%   The function saves the figure specified by the handle 'figHandle' to
%   disc. The function includes code that automatically removes obsolete
%   white margins that matlab puts around the figure by default. 
%
%
% Inputs:
%   figHandle: Handle to the figure containing the plots
%   fileName:  File name of the output file
%   options:   Optional arguments to the Matlab built-in print command.
%
% Outputs:
%
%
% Notes::
%   The code has been adopted from a resource that is published on 
%   Stackoverflow by the user Antonio. The link to the original code is:
%   http://stackoverflow.com/questions/3801730/get-rid-of-the-white-space-around-matlab-figures-pdf-output
%   Many thanks to Antonio for sharing this code!
%
% Examples::
%
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also createDataSheet, genericJoint, jointBuilder.      
%


% Match figure and paper units
set(figHandle,'PaperUnits','centimeters');

% Collect handles to children (e.g. subfigures)
a=get(figHandle,'Children');
nChild=length(a);

% Initialize variabels to manage individual children bounds.
% Bounds will contain lower-left and upper-right corners of plots plus one
% line to make sure single plots work
bounds=zeros(nChild+1,4);
bounds(end,1:2)=inf;
bounds(end,3:4)=-inf;

% Generate all coordinates of corners of graphs and store them in
% bounds as [lower-left-x lower-left-y upper-right-x upper-right-y] in
% the same unit system as paper (centimeters here)
for i=1:nChild
    set(a(i),'Unit','centimeters');
    pos=get(a(i),'Position');
    inset = [0 0 0 0] ;
    if isfield(a(i),'TightInset')
        inset=get(a(i),'TightInset');
    end
    bounds(i,:)=[pos(1)-inset(1) pos(2)-inset(2) ...
        pos(1)+pos(3)+inset(3) pos(2)+pos(4)+inset(4)];
end

% Compute the rectangular convex hull of all plots and store that info
% in mypos as [lower-left-x lower-left-y width height] in centimeters
auxmin=min(bounds(:,1:2));
auxmax=max(bounds(:,3:4));
mypos=[auxmin auxmax-auxmin];

% Set the paper to the exact size of the on-screen figure using
% figure property PaperSize [width height]
set(figHandle,'PaperSize',[mypos(3) mypos(4)]);

% Ensure that paper position mode is in manual in order for the
% printer driver to honor the figure properties
set(figHandle,'PaperPositionMode', 'manual');

% Use the PaperPosition four-element vector [left, bottom, width, height]
% to control the location on printed page; place it using horizontal and
% vertical negative offsets equal to the lower-left coordinates of the
% rectangular convex hull of the plot, and increase the size of the figure
% accordingly
set(figHandle,'PaperPosition',[-mypos(1) -mypos(2) ...
    mypos(3)+mypos(1) mypos(4)+mypos(2)]);

% Print stuff and pass optional parameters of the print command
if verLessThan('matlab','9.1') % check for backwards compatibility
    print('-dpdf',outfilename,varargin{:});
else
    print('-bestfit','-dpdf',outfilename,varargin{:});
end