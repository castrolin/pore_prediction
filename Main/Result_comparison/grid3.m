%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Grid Lines for 4-D, 5D or 6-D Scatter Plots     %
%              with MATLAB Implementation              %
%                                                      %
% Author: Ph.D. Eng. Hristo Zhivomirov        11/25/19 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grid3(varargin)
% function: grid3(vararg)
% varargin - type 'on' in the place of varargin if one want to add grid
%            lines to the current axes;
%          - type 'off' in the place of varargin if one want to remove the 
%            grid lines from the current axes.
% check for relevance
if isempty(findobj(gca, 'Type', 'Scatter'))
    error('The grid3 command operates only with 3D scatter plots.')
end
% check the varargin
if strcmp(varargin, 'on')
    r = 'r'; g = 'g'; b = 'b';
elseif strcmp(varargin, 'off')
    r = 'w'; g = 'w'; b = 'w';
else
    return
end
    
% remove all grid lines from the current axes
grid off
% set the axis limits to the range of the data
axis tight
% set the hold state to on
hold on
% obtain the Cartesian coordinate data
xdata = unique(get(findobj(gca, 'Type', 'Scatter'), 'XData'));
ydata = unique(get(findobj(gca, 'Type', 'Scatter'), 'YData'));
zdata = unique(get(findobj(gca, 'Type', 'Scatter'), 'ZData'));
minxdata = min(xdata);
maxxdata = max(xdata);
minydata = min(ydata);
maxydata = max(ydata);
minzdata = min(zdata); 
maxzdata = max(zdata);
% form x-grid and y-grid
for z = 1:length(zdata)
    % x-grid
    for x = 1:length(xdata)
        plot3([xdata(x) xdata(x)], ...
              [minydata maxydata], ...
              [zdata(z) zdata(z)], ':', 'Color', r, 'LineWidth', 0.125)
    end
    
    % y-grid
    for y = 1:length(ydata)
        plot3([minxdata maxxdata], ...
              [ydata(y) ydata(y)], ...
              [zdata(z) zdata(z)], ':', 'Color', g, 'LineWidth', 0.125)
    end
end
% form z-grid
for y = 1:length(ydata)
    for x = 1:length(xdata)
        plot3([xdata(x) xdata(x)], ...
              [ydata(y) ydata(y)], ...
              [minzdata maxzdata], ':', 'Color', b, 'LineWidth', 0.125)
    end
end
% set the hold state to off
hold off